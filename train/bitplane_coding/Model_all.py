import torch
import numpy as np

from Model_train_predict import Model_train_predict
from Model_train_update import Model_train_update
from Model_space_entropy import Model_space_entropy
#
alpha_level = np.sqrt(23.0/32.0)
beta_level = np.sqrt(3.0 / 2.0)
# alpha_level = 1.0
# beta_level = 1.0
class Model_all(torch.nn.Module):
    def __init__(self, wacelet_trainable_set, train_step):
        super(Model_all, self).__init__()

        self.preditct_wavelet_encode_H1 = Model_train_predict()
        self.update_wavelet_encode_H1 = Model_train_update()
        self.space_entropy_H1 = Model_space_entropy(wacelet_trainable_set, train_step)
        self.preditct_wavelet_decode_H1 = Model_train_predict()
        self.update_wavelet_decode_H1 = Model_train_update()  

        self.preditct_wavelet_encode_H2 = Model_train_predict()
        self.update_wavelet_encode_H2 = Model_train_update()
        self.space_entropy_H2 = Model_space_entropy(wacelet_trainable_set, train_step)
        self.preditct_wavelet_decode_H2 = Model_train_predict()
        self.update_wavelet_decode_H2 = Model_train_update()

        self.preditct_wavelet_encode_H3 = Model_train_predict()
        self.update_wavelet_encode_H3 = Model_train_update()
        self.space_entropy_H3 = Model_space_entropy(wacelet_trainable_set, train_step)
        self.preditct_wavelet_decode_H3 = Model_train_predict()  
        self.update_wavelet_decode_H3 = Model_train_update()

        self.space_entropy_L = Model_space_entropy(wacelet_trainable_set, train_step)
        self.part_bitplane = 6

    def forward(self, frame, mask, mv, grid_up, grid_org, train, alpha_list):
        Batch_size = frame.size()[0]
        frame = frame.unsqueeze(2)
        mask = mask.unsqueeze(2)
        

        H1_list = []
        for i in range(4):
            predict_john, _, predict_cnn = self.preditct_wavelet_encode_H1(frame[:, i * 2, :, :],
                                                                        mv[:, 4 * i:4 * i + 2, :, :],
                                                                        mask[:, i * 2, :, :],
                                                                        frame[:, i * 2 + 2, :, :],
                                                                        mv[:, 4 * i + 2:4 * i + 4, :, :],
                                                                        mask[:, i * 2 + 1, :, :], grid_up,
                                                                        grid_org)
            H_tmp = frame[:, i * 2 + 1, :, :] - predict_cnn
            H1_list.append(H_tmp*alpha_level)
        with torch.no_grad():
            for i in range(4, 4 + 3+ 4):
                predict_john, _, predict_cnn = self.preditct_wavelet_encode_H1(frame[:, i * 2, :, :],
                                                                               mv[:, 4 * i:4 * i + 2, :, :],
                                                                               mask[:, i * 2, :, :],
                                                                               frame[:, i * 2 + 2, :, :],
                                                                               mv[:, 4 * i + 2:4 * i + 4, :, :],
                                                                               mask[:, i * 2 + 1, :, :], grid_up,
                                                                               grid_org)
                H_tmp = frame[:, i * 2 + 1, :, :] - predict_cnn
                H1_list.append(H_tmp*alpha_level)
        H1 = torch.cat(H1_list[:4], 0)
        rec_H1, bits_H1, wavelet_scale_H1 = self.space_entropy_H1(H1, self.part_bitplane, train, alpha_list[:, 0])  # 重建时域第二帧
        with torch.no_grad():
            rec_H1_nogrid, bits_H1_nogrid, wavelet_scale_H1_nogrid = self.space_entropy_H1(H1_list[4], self.part_bitplane, train, alpha_list[:,0])


        L1_list = []
        i = 0
        L_tmp = self.update_wavelet_encode_H1(H1_list[0],
                                            -mv[:, 44+4 * i:44+4 * i + 2, :, :],
                                            mask[:, 22+i * 2, :, :],
                                            H1_list[0],
                                            -mv[:, 44+4 * i + 2:44+4 * i + 4, :, :],
                                            mask[:, 22+i * 2 + 1, :, :], grid_up,
                                            grid_org)
        L1 = (frame[:, 0, :, :] + L_tmp/alpha_level) * beta_level
        L1_list.append(L1)
        for i in range(1, 4):
             L_tmp = self.update_wavelet_encode_H1(H1_list[i-1],
                                                    -mv[:, 44+4 * i:44+4 * i + 2, :, :],
                                                    mask[:, 22+i * 2, :, :],
                                                    H1_list[i],
                                                    -mv[:, 44+4 * i + 2:44+4 * i + 4, :, :],
                                                    mask[:, 22+i * 2 + 1, :, :], grid_up,
                                                    grid_org)
             L1 = (frame[:, 2 * i, :, :] + L_tmp/alpha_level) * beta_level
             L1_list.append(L1)
        with torch.no_grad():
            for i in range(4, 4+3+4):
                L_tmp = self.update_wavelet_encode_H1(H1_list[i - 1],
                                                      -mv[:, 44 + 4 * i:44+4 * i + 2, :, :],
                                                      mask[:, 22 + i * 2, :, :],
                                                      H1_list[i],
                                                      -mv[:, 44 + 4 * i + 2:44+4 * i + 4, :, :],
                                                      mask[:, 22 + i * 2 + 1, :, :], grid_up,
                                                      grid_org)
                L1 = (frame[:, 2 * i, :, :] + L_tmp/alpha_level) * beta_level
                L1_list.append(L1)

        H2_list = []
        for i in range(2):
            predict_john, _, predict_cnn = self.preditct_wavelet_encode_H2(L1_list[2*i],
                                                                        mv[:, 44*2 + 4 * i:44*2 + 4 * i + 2, :, :],
                                                                        mask[:, 22*2 + i * 2, :, :],
                                                                        L1_list[2*i+2],
                                                                        mv[:, 44*2 + 4 * i + 2:44*2 + 4 * i + 4, :, :],
                                                                        mask[:, 22*2 + i * 2 + 1, :, :], grid_up,
                                                                        grid_org)
            H_tmp = L1_list[2*i+1] - predict_cnn
            H2_list.append(H_tmp*alpha_level)
        with torch.no_grad():
            for i in range(2, 2+1+2):
                predict_john, _, predict_cnn = self.preditct_wavelet_encode_H2(L1_list[2 * i],
                                                                               mv[:, 44 * 2 + 4 * i:44 * 2 + 4 * i + 2,
                                                                               :, :],
                                                                               mask[:, 22 * 2 + i * 2, :, :],
                                                                               L1_list[2 * i + 2],
                                                                               mv[:,
                                                                               44 * 2 + 4 * i + 2:44 * 2 + 4 * i + 4, :,
                                                                               :],
                                                                               mask[:, 22 * 2 + i * 2 + 1, :, :],
                                                                               grid_up,
                                                                               grid_org)
                H_tmp = L1_list[2 * i + 1] - predict_cnn
                H2_list.append(H_tmp*alpha_level)
        H2 = torch.cat(H2_list[:2], 0)
        rec_H2, bits_H2, wavelet_scale_H2 = self.space_entropy_H2(H2, self.part_bitplane, train, alpha_list[:, 1])  # 重建时域第二帧
        with torch.no_grad():
            rec_H2_nogrid, bits_H2_nogrid, wavelet_scale_H2_nogrid = self.space_entropy_H2(H2_list[2], self.part_bitplane, train, alpha_list[:,1])
                     

        L2_list = []
        i = 0
        L_tmp = self.update_wavelet_encode_H2(H2_list[0],
                                            -mv[:, 88+20+4 * i:88+20+4 * i + 2, :, :],
                                            mask[:, 44+10+i * 2, :, :],
                                            H2_list[0],
                                            -mv[:, 88+20+4 * i + 2:88+20+4 * i + 4, :, :],
                                            mask[:, 44+10+i * 2 + 1, :, :], grid_up,
                                            grid_org)
        L2 = (L1_list[2*i] + L_tmp/alpha_level)*beta_level
        L2_list.append(L2)
        i = 1
        L_tmp = self.update_wavelet_encode_H2(H2_list[i-1],
                                              -mv[:, 88 + 20 + 4 * i:88 + 20 + 4 * i + 2, :, :],
                                              mask[:, 44 + 10 + i * 2, :, :],
                                              H2_list[i],
                                              -mv[:, 88 + 20 + 4 * i + 2:88 + 20 + 4 * i + 4, :, :],
                                              mask[:, 44 + 10 + i * 2 + 1, :, :], grid_up,
                                              grid_org)
        L2 = (L1_list[2*i] + L_tmp/alpha_level)*beta_level
        L2_list.append(L2)
        with torch.no_grad():
            for i in range(2, 2+1+2):
                L_tmp = self.update_wavelet_encode_H2(H2_list[i - 1],
                                                      -mv[:, 88 + 20 + 4 * i:88 + 20 + 4 * i + 2, :, :],
                                                      mask[:, 44 + 10 + i * 2, :, :],
                                                      H2_list[i],
                                                      -mv[:, 88 + 20 + 4 * i + 2:88 + 20 + 4 * i + 4, :, :],
                                                      mask[:, 44 + 10 + i * 2 + 1, :, :], grid_up,
                                                      grid_org)
                L2 = (L1_list[2*i] + L_tmp/alpha_level)*beta_level
                L2_list.append(L2)

        H3_list = []
        i = 0
        predict_john, _, predict_cnn = self.preditct_wavelet_encode_H3(L2_list[2*i],
                                                                    mv[:, 128 + 4 * i:128 + 4 * i + 2, :, :],
                                                                    mask[:, 64 + 2 * i, :, :],
                                                                    L2_list[2*i+2],
                                                                    mv[:, 128 + 4 * i + 2:128 + 4 * i + 4, :, :],
                                                                    mask[:, 64 + 2 * i + 1, :, :], grid_up,
                                                                    grid_org)
        H_tmp = L2_list[2*i+1] - predict_cnn
        H3_list.append(H_tmp*alpha_level)
        with torch.no_grad():
            i = 1
            predict_john, _, predict_cnn = self.preditct_wavelet_encode_H3(L2_list[2 * i],
                                                                           mv[:, 128 + 4 * i:128 + 4 * i + 2, :, :],
                                                                           mask[:, 64 + 2 * i, :, :],
                                                                           L2_list[2 * i + 2],
                                                                           mv[:, 128 + 4 * i + 2:128 + 4 * i + 4, :, :],
                                                                           mask[:, 64 + 2 * i + 1, :, :], grid_up,
                                                                           grid_org)
            H_tmp = L2_list[2 * i + 1] - predict_cnn
            H3_list.append(H_tmp*alpha_level)
        H3_rec_list = []
        rec_H3, bits_H3, wavelet_scale_H3 = self.space_entropy_H3(H3_list[0], self.part_bitplane, train, alpha_list[:,2])  # 重建时域第二帧
        H3_rec_list.append(rec_H3)
        with torch.no_grad():
            rec_H3_nogrid, bits_H3_nogrid, wavelet_scale_H3_nogrid = self.space_entropy_H3(H3_list[1], self.part_bitplane, train, alpha_list[:,2])
            H3_rec_list.append(rec_H3_nogrid)

        L3_list = []
        i = 0
        L_tmp = self.update_wavelet_encode_H3(H3_list[0],
                                            -mv[:, 128+8+4 * i:128+8+4 * i + 2, :, :],
                                            mask[:, 64+4+i * 2, :, :],
                                            H3_list[0],
                                            -mv[:, 128+8+4 * i + 2:128+8+4 * i + 4, :, :],
                                            mask[:, 64+4+i * 2 + 1, :, :], grid_up,
                                            grid_org)
        L3 = (L2_list[0]+L_tmp/alpha_level)*beta_level
        L3_list.append(L3)
        with torch.no_grad():
            i = 1
            L_tmp = self.update_wavelet_encode_H3(H3_list[0],
                                                  -mv[:, 128 + 8 + 4 * i:128 + 8 + 4 * i + 2, :, :],
                                                  mask[:, 64 + 4 + i * 2, :, :],
                                                  H3_list[1],
                                                  -mv[:, 128 + 8 + 4 * i + 2:128 + 8 + 4 * i + 4, :, :],
                                                  mask[:, 64 + 4 + i * 2 + 1, :, :], grid_up,
                                                  grid_org)
            L3 = (L2_list[2] + L_tmp/alpha_level)*beta_level
            L3_list.append(L3)

        L3_rec_list = []
        rec_L, bits_L, wavelet_scale_L = self.space_entropy_L(L3_list[0], self.part_bitplane, train, alpha_list[:,3])  # 重建时域第1帧
        L3_rec_list.append(rec_L)
        with torch.no_grad():
            rec_L_nogrid, bits_L_nogrid, wavelet_scale_L_nogrid = self.space_entropy_L(L3_list[1], self.part_bitplane, train, alpha_list[:,3])  # 重建时域第1帧
            # rec_L_nogrid = self.image_enhance_L(rec_L_nogrid)
            L3_rec_list.append(rec_L_nogrid)

        L3_rec_update_list = []
        i = 0
        L_tmp = self.update_wavelet_decode_H3(H3_rec_list[0],
                                              -mv[:, 128 + 8 + 4 * i:128 + 8 + 4 * i + 2, :, :],
                                              mask[:, 64 + 4 + i * 2, :, :],
                                              H3_rec_list[0],
                                              -mv[:, 128 + 8 + 4 * i + 2:128 + 8 + 4 * i + 4, :, :],
                                              mask[:, 64 + 4 + i * 2 + 1, :, :], grid_up,
                                              grid_org)
        L3 = L3_rec_list[0]/beta_level - L_tmp/alpha_level
        L3_rec_update_list.append(L3)
        with torch.no_grad():
            i = 1
            L_tmp = self.update_wavelet_decode_H3(H3_rec_list[0],
                                                  -mv[:, 128 + 8 + 4 * i:128 + 8 + 4 * i + 2, :, :],
                                                  mask[:, 64 + 4 + i * 2, :, :],
                                                  H3_rec_list[1],
                                                  -mv[:, 128 + 8 + 4 * i + 2:128 + 8 + 4 * i + 4, :, :],
                                                  mask[:, 64 + 4 + i * 2 + 1, :, :], grid_up,
                                                  grid_org)
            L3 = L3_rec_list[1]/beta_level - L_tmp/alpha_level
            L3_rec_update_list.append(L3)

        i = 0
        predict_john, _, predict_cnn = self.preditct_wavelet_decode_H3(L3_rec_update_list[0],
                                                                       mv[:, 128 + 4 * i:128 + 4 * i + 2, :, :],
                                                                       mask[:, 64 + 2 * i, :, :],
                                                                       L3_rec_update_list[1],
                                                                       mv[:, 128 + 4 * i + 2:128 + 4 * i + 4, :, :],
                                                                       mask[:, 64 + 2 * i + 1, :, :],
                                                                       grid_up,
                                                                       grid_org)
        H_tmp = H3_rec_list[0]/alpha_level + predict_cnn
        # H_tmp = self.image_enhance_H3(H_tmp, predict_cnn, H3_rec_list[0]/alpha_level)
        L3_rec_update_list.insert(1, H_tmp)

        L2_rec_update_list = []
        i = 0
        L_tmp = self.update_wavelet_decode_H2(rec_H2[Batch_size * i:Batch_size * (i + 1)],
                                              -mv[:, 88 + 20 + 4 * i:88 + 20 + 4 * i + 2, :, :],
                                              mask[:, 44 + 10 + i * 2, :, :],
                                              rec_H2[Batch_size * i:Batch_size * (i + 1)],
                                              -mv[:, 88 + 20 + 4 * i + 2:88 + 20 + 4 * i + 4, :, :],
                                              mask[:, 44 + 10 + i * 2 + 1, :, :], grid_up,
                                              grid_org)
        L2 = L3_rec_update_list[i]/beta_level - L_tmp/alpha_level
        L2_rec_update_list.append(L2)
        i = 1
        L_tmp = self.update_wavelet_decode_H2(rec_H2[Batch_size * (i-1):Batch_size * i],
                                              -mv[:, 88 + 20 + 4 * i:88 + 20 + 4 * i + 2, :, :],
                                              mask[:, 44 + 10 + i * 2, :, :],
                                              rec_H2[Batch_size * i:Batch_size * (i + 1)],
                                              -mv[:, 88 + 20 + 4 * i + 2:88 + 20 + 4 * i + 4, :, :],
                                              mask[:, 44 + 10 + i * 2 + 1, :, :], grid_up,
                                              grid_org)
        L2 = L3_rec_update_list[i]/beta_level - L_tmp/alpha_level
        L2_rec_update_list.append(L2)
        with torch.no_grad():
            i = 2
            L_tmp = self.update_wavelet_decode_H2(rec_H2[Batch_size * (i-1):Batch_size * i],
                                                  -mv[:, 88 + 20 + 4 * i:88 + 20 + 4 * i + 2, :, :],
                                                  mask[:, 44 + 10 + i * 2, :, :],
                                                  rec_H2_nogrid,
                                                  -mv[:, 88 + 20 + 4 * i + 2:88 + 20 + 4 * i + 4, :, :],
                                                  mask[:, 44 + 10 + i * 2 + 1, :, :], grid_up,
                                                  grid_org)
            L2 = L3_rec_update_list[i]/beta_level - L_tmp/alpha_level
            L2_rec_update_list.append(L2)

        for i in range(2):
            predict_john, _, predict_cnn = self.preditct_wavelet_decode_H2(L2_rec_update_list[i+i],
                                                                           mv[:, 44 * 2 + 4 * i:44 * 2 + 4 * i + 2, :,
                                                                           :],
                                                                           mask[:, 22 * 2 + i * 2, :, :],
                                                                           L2_rec_update_list[i + 1+i],
                                                                           mv[:, 44 * 2 + 4 * i + 2:44 * 2 + 4 * i + 4,
                                                                           :, :],
                                                                           mask[:, 22 * 2 + i * 2 + 1, :, :], grid_up,
                                                                           grid_org)
            H_tmp = rec_H2[Batch_size * i:Batch_size * (i + 1)]/alpha_level + predict_cnn
            # H_tmp = self.image_enhance_H2(H_tmp, predict_cnn, rec_H2[Batch_size * i:Batch_size * (i + 1)]/alpha_level)
            L2_rec_update_list.insert(i+1+i, H_tmp)


        L1_rec_update_list = []
        i = 0
        L_tmp = self.update_wavelet_decode_H1(rec_H1[Batch_size * i:Batch_size * (i + 1)],
                                            -mv[:, 44+4 * i:44+4 * i + 2, :, :],
                                            mask[:, 22+i * 2, :, :],
                                            rec_H1[Batch_size * i:Batch_size * (i + 1)],
                                            -mv[:, 44+4 * i + 2:44+4 * i + 4, :, :],
                                            mask[:, 22+i * 2 + 1, :, :], grid_up,
                                            grid_org)
        L1 = L2_rec_update_list[i]/beta_level - L_tmp/alpha_level
        L1_rec_update_list.append(L1)
        for i in range(1, 4):
            L_tmp = self.update_wavelet_decode_H1(rec_H1[Batch_size * (i-1):Batch_size * i],
                                                -mv[:, 44+4 * i:44+4 * i + 2, :, :],
                                                mask[:, 22+i * 2, :, :],
                                                rec_H1[Batch_size * i:Batch_size * (i+1)],
                                                -mv[:, 44+4 * i + 2:44+4 * i + 4, :, :],
                                                mask[:, 22+i * 2 + 1, :, :], grid_up,
                                                grid_org)
            L1 = L2_rec_update_list[i]/beta_level - L_tmp/alpha_level
            L1_rec_update_list.append(L1)
        with torch.no_grad():
            i = 4
            L_tmp = self.update_wavelet_decode_H1(rec_H1[Batch_size * (i-1):Batch_size * i],
                                                  -mv[:, 44 + 4 * i:44 + 4 * i + 2, :, :],
                                                  mask[:, 22 + i * 2, :, :],
                                                  rec_H1_nogrid,
                                                  -mv[:, 44 + 4 * i + 2:44 + 4 * i + 4, :, :],
                                                  mask[:, 22 + i * 2 + 1, :, :], grid_up,
                                                  grid_org)
            L1 = L2_rec_update_list[i]/beta_level - L_tmp/alpha_level
            L1_rec_update_list.append(L1)

        for i in range(4):
            predict_john, _, predict_cnn = self.preditct_wavelet_decode_H1(L1_rec_update_list[i+i],
                                                                        mv[:, 4 * i:4 * i + 2, :, :],
                                                                        mask[:, i * 2, :, :],
                                                                        L1_rec_update_list[i+1+i],
                                                                        mv[:, 4 * i + 2:4 * i + 4, :, :],
                                                                        mask[:, i * 2 + 1, :, :], grid_up,
                                                                        grid_org)
            H_tmp = rec_H1[Batch_size * i:Batch_size * (i + 1)]/alpha_level + predict_cnn
            # H_tmp = self.image_enhance_H1(H_tmp, predict_cnn, rec_H1[Batch_size * i:Batch_size * (i + 1)]/alpha_level)
            L1_rec_update_list.insert(i + 1 + i, H_tmp)
        return torch.cat(L1_rec_update_list,1), [bits_H1, bits_H2, bits_H3, bits_L], [wavelet_scale_H1,
               wavelet_scale_H2, wavelet_scale_H3, wavelet_scale_L]