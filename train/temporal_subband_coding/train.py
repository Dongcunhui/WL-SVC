import torch
import argparse
from torch.autograd import Variable
import os
from torch.utils.data import DataLoader
import torch.optim as optim
import datetime
import random
import numpy as np

from Model_all import Model_all
from Mydata import MyData

# os.environ['CUDA_VISIBLE_DEVICES']='2'
os.environ["NCCL_DEBUG"] = "INFO"
parser = argparse.ArgumentParser(description='john-end2end')
parser.add_argument('--input_train_dir', type=str, default=r'./name_list.txt')
parser.add_argument('--input_test_dir', type=str, default=r'./test/')
parser.add_argument('--checkpoint_dir', type=str, default='./model_dir/')
parser.add_argument('--log_dir', type=str, default='./')
parser.add_argument('--epochs', type=int, default=1000)
parser.add_argument('--batch_size', type=int, default=3)
parser.add_argument('--load_model', type=str,default=None)
parser.add_argument('--rate_b', type=str,default=None)



def setup_seed(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True

def to_variable(x):
    if torch.cuda.is_available():
        x = x.cuda()
    return Variable(x)

def main():
    args = parser.parse_args()
    input_train_dir = args.input_train_dir
    input_test_dir = args.input_test_dir
    ckpt_dir = args.checkpoint_dir
    if not os.path.exists(ckpt_dir):
        os.makedirs(ckpt_dir)
    if not os.path.exists(args.log_dir):
        os.makedirs(args.log_dir)
    logfile = open(args.log_dir + '/log.txt', 'w')
    total_epoch = args.epochs
    batch_size = args.batch_size

    wacelet_trainable_set = False
    lambda_d = 0.064
    temporal_number = 3
    frame_number = 2 ** temporal_number
    alpha_list = [10, 10, 10, 10, 2]  # not useful
    crop_size = 256
    train_step = 4
    model = Model_all(wacelet_trainable_set, train_step)
    epoch = 1
    accumulation_steps = 1
    Optimizer = optim.Adam(filter(lambda p: p.requires_grad, model.parameters()), lr=1*1e-4)

    logfile.flush()
    dataset_train = MyData(input_train_dir, crop_size, train=1)
    dataset_test = MyData(input_test_dir, crop_size, train=0)

    if args.load_model is not None:
        checkpoint = torch.load(args.load_model)
        state_dict = checkpoint['state_dict']
        model.load_state_dict(state_dict)
        epoch = checkpoint['epoch']
        alpha_list = checkpoint['alpha_list']
        Optimizer.load_state_dict(checkpoint['Optimizer'])
        if torch.cuda.is_available():
            for state in Optimizer.state.values():
                for k, v in state.items():
                    if torch.is_tensor(v):
                        state[k] = v.cuda()
        print('Load pre-trained model [' + args.load_model + '] succeed!')
        logfile.write('Load pre-trained model [' + args.load_model + '] succeed!' + '\n')
        logfile.flush()
    else:
        if args.rate_b is not None:
            model_dict = model.state_dict()
            checkpoint = torch.load(args.rate_b)
            part_dict_org = checkpoint['state_dict']
            print("pre train args.rate_b", end=" ")
            logfile.write("pre train args.rate_b" + ' ')
            part_dict = {k: v for k, v in part_dict_org.items() if "space_entropy" in k and "_E" not in k and "scale" not in k}
            for name, param in part_dict.items():
                print("pretrain "+name, end=" ")
                logfile.write("pretrain "+name + ' ')
            model_dict.update(part_dict)
            model.load_state_dict(model_dict)
            epoch = 1
            print('Load pre-trained part model [' +args.rate_b + '] succeed!')
            logfile.write('Load pre-trained part model [' + args.rate_b + '] succeed!' + '\n')
            logfile.flush()

        
    alpha_list = np.array(alpha_list, dtype=np.float32)
    alpha_list = torch.from_numpy(alpha_list)
    if torch.cuda.is_available():
        alpha_list = alpha_list.cuda()
        alpha_list = alpha_list.unsqueeze(0)
        alpha_list = torch.repeat_interleave(alpha_list, repeats=torch.cuda.device_count(), dim=0)

    if torch.cuda.is_available():
        if torch.cuda.device_count() > 1:
            print("Let's use", torch.cuda.device_count(), "GPUs!")
            model = torch.nn.DataParallel(model).cuda()
        else:
            model = model.cuda()
            
    train_loader = DataLoader(dataset=dataset_train, batch_size=batch_size, shuffle=True, num_workers=8, pin_memory=True)
    test_loader = DataLoader(dataset=dataset_test, batch_size=1, shuffle=True, num_workers=2, pin_memory=True)
    max_step_train = train_loader.__len__()
    max_step_test = test_loader.__len__()
    
    print("can train parameters:", end=' ')
    logfile.write("can train parameters:" + ' ')
    for name, param in model.named_parameters():
        if param.requires_grad:
            print(name, end=' ')
            logfile.write(name + ' ')
    print(' ')
    logfile.write('\n')
    logfile.flush()

    scheduler = torch.optim.lr_scheduler.StepLR(Optimizer, step_size=1, gamma=0.8)
    mse_func = torch.nn.MSELoss()
    
    while( epoch < total_epoch):
        Optimizer.zero_grad()
        print(str([alpha.item() for alpha in alpha_list[0]]))
        logfile.write('alpha_L: ' + str([alpha.item() for alpha in alpha_list[0]]) + '\n')
        logfile.flush()
        model.train()

        psnr_list_all = []
        bpp_list_all = []
        train_scale_list_all =[]
        for i in range(frame_number):
            psnr_list_all.append(0.0)
        for i in range(10):
            bpp_list_all.append(0.0)
            train_scale_list_all.append(0.0)
        loss_all = 0.
        psnr_all = 0.
        bpp_all = 0.

        W = crop_size * 4
        H = crop_size * 4
        xx = torch.arange(0, W).view(1, -1).repeat(H, 1)
        yy = torch.arange(0, H).view(-1, 1).repeat(1, W)
        xx = xx.view(1, 1, H, W).repeat(batch_size, 1, 1, 1)
        yy = yy.view(1, 1, H, W).repeat(batch_size, 1, 1, 1)
        grid_up = torch.cat((xx, yy), 1).float()
        grid_up = to_variable(grid_up)

        W = crop_size
        H = crop_size
        xx = torch.arange(0, W).view(1, -1).repeat(H, 1)
        yy = torch.arange(0, H).view(-1, 1).repeat(1, W)
        xx = xx.view(1, 1, H, W).repeat(batch_size, 1, 1, 1)
        yy = yy.view(1, 1, H, W).repeat(batch_size, 1, 1, 1)
        grid_org = torch.cat((xx, yy), 1).float()
        grid_org = to_variable(grid_org)

        for batch_idx, input_img in enumerate(train_loader):
            input_img_v = to_variable(input_img)
            size = input_img_v.size()
            B = size[0]
            grid_org = grid_org[:B]
            grid_up = grid_up[:B]
            frame = input_img_v[:, :23, :, :]
            mask = input_img_v[:, 23:95, :, :]
            mv = input_img_v[:, 95:, :, :] / 16.0
            # Optimizer.zero_grad()
            rec, bits_list, train_scale_list = model(frame, mask, mv, grid_up, grid_org, 1, alpha_list)
            bits = torch.sum(sum(bits_list))
            bpp = bits / size[0] / size[2] / size[3] / frame_number
            psnr_list = []
            for frame_order in range(frame_number):
                mse = mse_func(rec[:,frame_order,:,:], frame[:,frame_order,:,:])
                psnr = 10. * torch.log10(255. * 255. / mse)
                psnr_list.append(psnr)
            psnr_ave = torch.sum(sum(psnr_list))/frame_number
            loss = (bpp)/accumulation_steps

            loss.backward()
            if ((batch_idx+1) % accumulation_steps) == 0:
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=20, norm_type=2)
                Optimizer.step()
                Optimizer.zero_grad()

            for i in range(frame_number):
                psnr_list_all[i] += psnr_list[i].item()
            for i in range(len(bits_list)):
                bpp_list_all[i] += (torch.sum(bits_list[i]) / size[0] / size[2] / size[3] / frame_number).item()
            for i in range(len(train_scale_list)):    
                train_scale_list_all[i] += (torch.sum(train_scale_list[i]) / torch.cuda.device_count()).item()
            psnr_all = psnr_all + psnr_ave.item()
            bpp_all = bpp_all + bpp.item()
            loss_all = loss_all + loss.item() * accumulation_steps

            if batch_idx % 100 == 0:
                logfile.write('Train Epoch: [' + str(epoch) + '/' + str(total_epoch) + ']   ' + 'Step: [' + str(
                    batch_idx) + '/' + str(max_step_train) + ']   '  + 'train loss: ' + str(loss.item() * accumulation_steps) + '/' + '/' + str(
                    bpp.item()) + '/' + str(psnr_ave.item()) + '\n')
                logfile.flush()
                print('Train Epoch: [' + str(epoch) + '/' + str(total_epoch) + ']   ' + 'Step: [' + str(
                    batch_idx) + '/' + str(max_step_train) + ']   ' + 'train loss: ' + str(loss.item() * accumulation_steps) + '/' + '/' + str(
                    bpp.item()) + '/' + str(psnr_ave.item()) + '\n')

        for i in range(frame_number):
            psnr_list_all[i] = psnr_list_all[i] / max_step_train
        for i in range(len(train_scale_list_all)):
            train_scale_list_all[i] = train_scale_list_all[i] / max_step_train
        for i in range(len(bpp_list_all)):    
            bpp_list_all[i] = bpp_list_all[i] / max_step_train
        logfile.write('true train loss_mean: ' + str(loss_all / max_step_train) + '\n')
        logfile.write('true train bpp_mean: ' + str(bpp_all / max_step_train) + '\n')
        logfile.write('true train psnr_mean: ' + str(psnr_all / max_step_train) + '\n')
        logfile.write('true train psnr_list: ' + str(psnr_list_all) + '\n')
        logfile.write('true train bpp_list: ' + str(bpp_list_all) + '\n')
        logfile.write('true train train_scale: ' + str(train_scale_list_all) + '\n')
        logfile.flush()
        print('true train loss_mean: ' + str(loss_all / max_step_train))
        print('true train bpp_mean: ' + str(bpp_all / max_step_train))
        print('true train psnr_mean: ' + str(psnr_all / max_step_train))
        print('true train psnr_list: ' + str(psnr_list_all))
        print('true train bpp_list: ' + str(bpp_list_all))
        print('true train train_scale: ' + str(train_scale_list_all))

        if epoch % 1 == 0:
            if torch.cuda.device_count() > 1:
                torch.save({'epoch': epoch, 'alpha_list': [alpha.item() for alpha in alpha_list[0]], 'state_dict': model.module.state_dict(), 'Optimizer' : Optimizer.state_dict()},
                       ckpt_dir + '/epoch' +  str(epoch//30*30).zfill(3) + '.pth' ,_use_new_zipfile_serialization=False)
            else:
                torch.save({'epoch': epoch, 'alpha_list': [alpha.item() for alpha in alpha_list[0]], 'state_dict': model.state_dict(),'Optimizer' : Optimizer.state_dict()},
                       ckpt_dir + '/epoch' +  str(epoch//30*30).zfill(3) + '.pth',_use_new_zipfile_serialization=False)

        if epoch % 5 == 0:
            # test
            psnr_list_all = []
            bpp_list_all = []
            train_scale_list_all = []
            for i in range(frame_number):
                psnr_list_all.append(0.0)
            for i in range(10):
                bpp_list_all.append(0.0)
                train_scale_list_all.append(0.0)
            loss_all = 0.
            psnr_all = 0.
            bpp_all = 0.
            model.eval()

            for batch_idx, input_img in enumerate(test_loader):
                with torch.no_grad():
                    size = input_img.size()
                    width = size[3]
                    height = size[2]
                    height = int(np.floor(height / 16)) * 16  # 16的整数倍
                    width = int(np.floor(width / 16)) * 16
                    input_img = input_img[:, :, :height, :width]

                    input_img_v = to_variable(input_img)
                    size = input_img_v.size()
                    width = size[3]
                    height = size[2]
                    W = width * 4
                    H = height * 4
                    xx = torch.arange(0, W).view(1, -1).repeat(H, 1)
                    yy = torch.arange(0, H).view(-1, 1).repeat(1, W)
                    xx = xx.view(1, 1, H, W).repeat(1, 1, 1, 1)
                    yy = yy.view(1, 1, H, W).repeat(1, 1, 1, 1)
                    grid_up = torch.cat((xx, yy), 1).float()
                    grid_up = to_variable(grid_up)

                    W = width
                    H = height
                    xx = torch.arange(0, W).view(1, -1).repeat(H, 1)
                    yy = torch.arange(0, H).view(-1, 1).repeat(1, W)
                    xx = xx.view(1, 1, H, W).repeat(1, 1, 1, 1)
                    yy = yy.view(1, 1, H, W).repeat(1, 1, 1, 1)
                    grid_org = torch.cat((xx, yy), 1).float()
                    grid_org = to_variable(grid_org)

                    frame = input_img_v[:, :23, :, :]
                    mask = input_img_v[:, 23:95, :, :]
                    mv = input_img_v[:, 95:, :, :] / 16.0

                    rec, bits_list, train_scale_list = model(frame, mask, mv, grid_up, grid_org, 0, alpha_list)
                    bits = torch.sum(sum(bits_list))
                    bpp = bits / size[0] / size[2] / size[3] / frame_number
                    psnr_list = []
                    for frame_order in range(frame_number):
                        mse = mse_func(rec[:, frame_order, :, :], frame[:, frame_order, :, :])
                        psnr = 10. * torch.log10(255. * 255. / mse)
                        psnr_list.append(psnr)
                    psnr_ave = torch.sum(sum(psnr_list)) / frame_number


                    loss = bpp

                    for i in range(frame_number):
                        psnr_list_all[i] += psnr_list[i].item()
                    for i in range(len(bits_list)):
                        bpp_list_all[i] += (torch.sum(bits_list[i]) / size[0] / size[2] / size[3] / frame_number).item()
                    for i in range(len(train_scale_list)):
                        train_scale_list_all[i] += (torch.sum(train_scale_list[i]) / torch.cuda.device_count()).item()
                    psnr_all = psnr_all + psnr_ave.item()
                    bpp_all = bpp_all + bpp.item()
                    loss_all = loss_all + loss.item()


            for i in range(frame_number):
                psnr_list_all[i] = psnr_list_all[i] / max_step_test
            for i in range(len(train_scale_list_all)):
                train_scale_list_all[i] = train_scale_list_all[i] / max_step_test
            for i in range(len(bpp_list_all)):
                bpp_list_all[i] = bpp_list_all[i] / max_step_test
            logfile.write('test test loss_mean: ' + str(loss_all / max_step_test) + '\n')
            logfile.write('test test bpp_mean: ' + str(bpp_all / max_step_test) + '\n')
            logfile.write('test test psnr_mean: ' + str(psnr_all / max_step_test) + '\n')
            logfile.write('test test psnr_list: ' + str(psnr_list_all) + '\n')
            logfile.write('test test bpp_list: ' + str(bpp_list_all) + '\n')
            logfile.write('test test train_scale: ' + str(train_scale_list_all) + '\n')
            logfile.flush()
            print('test test loss_mean: ' + str(loss_all / max_step_test))
            print('test test bpp_mean: ' + str(bpp_all / max_step_test))
            print('test test psnr_mean: ' + str(psnr_all / max_step_test))
            print('test test psnr_list: ' + str(psnr_list_all))
            print('test test bpp_list: ' + str(bpp_list_all))
            print('test test train_scale: ' + str(train_scale_list_all))

            if scheduler.get_lr()[0] > 1e-4:
                logfile.flush()
                scheduler.step()
            print("learn rate:", scheduler.get_lr()[0])
            logfile.write("learn rate:" + str(scheduler.get_lr()[0]) + '\n')
        print(datetime.datetime.now().strftime('%Y-%m-%d-%H:%M:%S'))
        logfile.write(str(datetime.datetime.now().strftime('%Y-%m-%d-%H:%M:%S')) + '\n')
        logfile.flush()
        if epoch % 20 == 0:
            for i in range(temporal_number):
                alpha_list[:,i] += 10.0
            alpha_list[:,temporal_number] += 1.0
        epoch = epoch + 1

    logfile.close()

if __name__ == "__main__":
    # setup_seed(20)
    main()

