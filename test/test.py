import torch
import argparse
from torch.autograd import Variable
import os
import datetime
from torch.nn import functional as F
import numpy as np

from Model_all_test import Model_all

os.environ['CUDA_VISIBLE_DEVICES'] = '1'
os.environ["NCCL_DEBUG"] = "INFO"
parser = argparse.ArgumentParser(description='WL-SVC')

def to_variable(x):
    if torch.cuda.is_available():
        x = x.cuda()
    return Variable(x)

def ndarray_nearest_neighbour_scaling(label, scale):
    assert len(label.shape) == 3
    c,h,w = label.shape

    label_new = np.zeros([c, h*scale, w*scale], dtype=label.dtype)


    y_pos = np.arange(h*scale)
    x_pos = np.arange(w*scale)
    y_pos = np.floor(y_pos / scale).astype(np.int32)
    x_pos = np.floor(x_pos / scale).astype(np.int32)

    y_pos = y_pos.reshape(y_pos.shape[0], 1)
    x_pos = x_pos.reshape(1, x_pos.shape[0])
    y_pos = np.tile(y_pos, (1, w*scale))
    x_pos = np.tile(x_pos, (h*scale, 1))
    assert y_pos.shape == x_pos.shape

    label_new[:, :, :] = label[:, y_pos[:, :], x_pos[:, :]]
    return label_new


def main():
    logfile = open('./log.txt', 'w')

    bitplane_number = 9
    train_step = 4
    model = Model_all(bitplane_number, train_step)
    logfile.flush()
    name_list = ["BasketballDrill_832x480_50"]
    path = r"data/"

    # load encode
    model_dict = model.state_dict()
    checkpoint = torch.load('model_all_encode.pth')
    part_dict = checkpoint['state_dict']
    part_dict = {k: v for k, v in part_dict.items() if "space_en" in k and "scale" not in k}
    model_dict.update(part_dict)
    model.load_state_dict(model_dict)
    discard_list = [0] # discard wavelet subband
    for name_sequ in name_list:
        for bitplane_encode_number in range(7, 2, -1): # encode bitplane number
            print("bitplane_encode_number:", bitplane_encode_number)
            for discard_number in discard_list:
                print("discard_number:", discard_number)
                # load post
                model_dict = model.state_dict()
                checkpoint = torch.load("wave_post/"+str(bitplane_encode_number-1)+'.pth')  # load learnable temporal wavelet inverse transform pre-train model
                part_dict = checkpoint['state_dict']
                part_dict = { k: v for k, v in part_dict.items() if "_wavelet_" in k}
                model_dict.update(part_dict)
                model.load_state_dict(model_dict)

                if torch.cuda.is_available():
                    model = model.cuda()

                mse_func = torch.nn.MSELoss()

                temporal_number = 3
                gop_size = 2**temporal_number
                channel_list = ["y", "u", "v"]
                bits_all_channel = []
                for i in range(1000):
                    bits_all_channel.append(0.0)
                rec_list = []
                for channel in channel_list:
                    print(name_sequ+"_"+channel)
                    model.eval()
                    with torch.no_grad():
                        frame = np.load(path+name_sequ+"_frame_"+channel+".npy")
                        mask = np.load(path + name_sequ + "_mask_"+channel+".npy")
                        mv = np.load(path + name_sequ + "_mv_"+channel+".npy") / 16.0  # mv is converted from float to int
                        frame = np.array(frame, dtype=np.float32)
                        frame = torch.from_numpy(frame).unsqueeze(0)
                        mask = np.array(mask, dtype=np.float32)
                        mask = torch.from_numpy(mask).unsqueeze(0)
                        mv = np.array(mv, dtype=np.float32)
                        mv = torch.from_numpy(mv).unsqueeze(0)
                        print("load date finish!")

                        size = frame.size()
                        width = size[3]
                        height = size[2]
                        frame_number = size[1]
                        gop_number = np.int(np.ceil(frame_number / gop_size))
                        pad = torch.zeros(frame.size()[0], (gop_number+1) * 8 - frame_number,
                                          frame.size()[2], frame.size()[3])
                        frame = torch.cat((frame, pad), 1)

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

                        rec, bits_list = model(frame, mask, mv, grid_up, grid_org, gop_number, bitplane_encode_number, discard_number)
                        rec_list.append(torch.clamp(torch.round(rec), 0, 255))
                        psnr_list = []
                        for frame_order in range(frame_number):
                            mse = mse_func(torch.clamp(torch.round(rec[:, frame_order, :, :]), 0, 255), frame[:, frame_order, :, :])
                            psnr = 10. * torch.log10(255. * 255. / mse)
                            psnr_list.append(psnr)
                            print(channel, frame_order, psnr.item())
                        psnr_ave = torch.sum(sum(psnr_list)) / frame_number
                        print(channel, "psnr_ave:", psnr_ave.item())
                        bits_all = 0
                        count = 0
                        for i in range(len(bits_list)):
                            for m in range(len(bits_list[i])):
                                bits_all += (bits_list[i][m].item())
                                print(i, m, bits_list[i][m].item())
                                bits_all_channel[count] += bits_list[i][m].item()
                                count += 1
                        print(channel, "bits_all:", bits_all/frame_number/width/height, bits_all)
                        # exit()
                all_channel_bits = 0
                for i in range(count):
                    print("channel_bits:", bits_all_channel[i]/8)
                    all_channel_bits += bits_all_channel[i]/8
                print("all_channel_bits:", all_channel_bits)

            # print(rec_list[0].size())
            # print("write rec yuv")
            # outpath = "rec_video/"
            # if not os.path.exists(outpath):
            #     os.makedirs(outpath)
            # fileout = open(outpath+name_sequ+"_"+str(rate)+".yuv", 'wb')
            # for frame_i in range(frame_number):
            #     tmp = rec_list[0].cpu().numpy()[0, frame_i, :, :]
            #     for height_i in range(tmp.shape[0]):
            #         for width_i in range(tmp.shape[1]):
            #             fileout.write(struct.pack('B', np.uint8(tmp[height_i, width_i])))
            #     tmp = rec_list[1].cpu().numpy()[0, frame_i, :, :]
            #     for height_i in range(tmp.shape[0]):
            #         for width_i in range(tmp.shape[1]):
            #             fileout.write(struct.pack('B', np.uint8(tmp[height_i, width_i])))
            #     tmp = rec_list[2].cpu().numpy()[0, frame_i, :, :]
            #     for height_i in range(tmp.shape[0]):
            #         for width_i in range(tmp.shape[1]):
            #             fileout.write(struct.pack('B', np.uint8(tmp[height_i, width_i])))
            # fileout.close()
        print(datetime.datetime.now().strftime('%Y-%m-%d-%H:%M:%S'))
        logfile.write(str(datetime.datetime.now().strftime('%Y-%m-%d-%H:%M:%S')) + '\n')
        logfile.flush()

    logfile.close()


if __name__ == "__main__":
    # setup_seed(20)
    main()