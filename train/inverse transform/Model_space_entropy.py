import torch
import torch.optim as optim
from torch.autograd import Variable
import math
import sys
from torch.nn import functional as F
import Quant            
import PixelCNN
import torch.nn as nn
import Wavelet_trans
import numpy as np
import downsample as downsample

class Model_space_entropy(torch.nn.Module):
    def __init__(self, wavelet_trainable_set, train_step):
        super(Model_space_entropy, self).__init__()

        self.trans_steps = train_step
        # self.coding_LL = PixelCNN.PixelCNN()
        # self.coding_HL_list = torch.nn.ModuleList([PixelCNN.PixelCNN_Context(1) for _i in range(self.trans_steps)])
        # self.coding_LH_list = torch.nn.ModuleList([PixelCNN.PixelCNN_Context(2) for _i in range(self.trans_steps)])
        # self.coding_HH_list = torch.nn.ModuleList([PixelCNN.PixelCNN_Context(3) for _i in range(self.trans_steps)])
        # for p in self.parameters():
            # p.requires_grad = False
            
        self.wavelet_transform = Wavelet_trans.Wavelet()
        self.quant = Quant.Quant()
        
        # self.coding_LL_E = PixelCNN.PixelCNN_Context(16+1)
        # self.coding_HL_list_E = torch.nn.ModuleList([PixelCNN.PixelCNN_Context(1+16+1) for _i in range(self.trans_steps)])
        # self.coding_LH_list_E = torch.nn.ModuleList([PixelCNN.PixelCNN_Context(2+16+1) for _i in range(self.trans_steps)])
        # self.coding_HH_list_E = torch.nn.ModuleList([PixelCNN.PixelCNN_Context(3+16+1) for _i in range(self.trans_steps)])

        # self.down1 = torch.nn.ModuleList([downsample.down_sample(1, 16) for _i in range(3)])
        # self.down2 = torch.nn.ModuleList([downsample.down_sample(2, 16) for _i in range(3)])
        # self.down3 = torch.nn.ModuleList([downsample.down_sample(3, 16) for _i in range(3)])
        # self.down4 = torch.nn.ModuleList([downsample.down_sample(4, 16) for _i in range(4)])
        # self.down = torch.nn.ModuleList([self.down1, self.down2, self.down3, self.down4])

    def forward(self, x, part_bitplane, train, alpha):
        # forward transform
        if not train:
            size = x.size()
            width = size[3]
            height = size[2]
            pad_h = int(np.ceil(height / 16)) * 16 - height  # 16的整数倍
            pad_w = int(np.ceil(width / 16)) * 16 - width
            paddings = (0, pad_w, 0, pad_h)
            x = F.pad(x, paddings, 'replicate')
        
        # 模拟前面一层bitplane
        self.scale = 2 ** (8 - part_bitplane + 1)
        LL = x
        HL_list_base = []
        LH_list_base = []
        HH_list_base = []
        for i in range(self.trans_steps):
            LL, HL, LH, HH = self.wavelet_transform.forward_trans(LL)
            if train:
                HL_list_base.append(self.quant(HL, self.scale) * self.scale)
                LH_list_base.append(self.quant(LH, self.scale) * self.scale)
                HH_list_base.append(self.quant(HH, self.scale) * self.scale)
            else:
                HL_list_base.append(self.quant(HL, self.scale) * self.scale)
                LH_list_base.append(self.quant(LH, self.scale) * self.scale)
                HH_list_base.append(self.quant(HH, self.scale) * self.scale)

        if train:
            LL_base = self.quant(LL, self.scale) * self.scale
        else:
            LL_base = self.quant(LL, self.scale) * self.scale
        LL = LL_base
        for i in range(self.trans_steps):
            j = self.trans_steps - 1 - i
            LL = self.wavelet_transform.inverse_trans(LL, HL_list_base[j], LH_list_base[j], HH_list_base[j])
        base_rec = LL
        
        # 当前bitplane
        self.scale = 2 ** (8 - part_bitplane)
        LL = x
        HL_list = []
        LH_list = []
        HH_list = []
        for i in range(self.trans_steps):
            LL, HL, LH, HH = self.wavelet_transform.forward_trans(LL)
            if train:
                HL_list.append(self.quant(HL, self.scale))
                LH_list.append(self.quant(LH, self.scale))
                HH_list.append(self.quant(HH, self.scale))
            else:
                HL_list.append(self.quant(HL, self.scale))
                LH_list.append(self.quant(LH, self.scale))
                HH_list.append(self.quant(HH, self.scale))

        if train:
            LL = self.quant(LL, self.scale)
        else:
            LL = self.quant(LL, self.scale)  
        
        # bits = self.coding_LL_E(LL - LL_base / self.scale, torch.cat((self.down[3][3](base_rec), LL_base/self.scale), 1))
        sign = torch.sign(LL)
        tmp = LL * sign + 0.5
        tmp[sign == 0] = 0
        LL = tmp * sign
        LL = (LL * self.scale)

        for i in range(self.trans_steps):
            j = self.trans_steps - 1 - i
            
            # bits = bits + self.coding_HL_list_E[j](HL_list[j]-HL_list_base[j]/self.scale,torch.cat((LL/self.scale, HL_list_base[j]/self.scale, self.down[j][0](base_rec)),1))
            sign = torch.sign(HL_list[j])
            tmp = HL_list[j] * sign + 0.5
            tmp[sign == 0] = 0
            HL_list[j] = tmp * sign
            HL_list[j] = (HL_list[j] * self.scale)
            
            # bits = bits + self.coding_LH_list_E[j](LH_list[j]-LH_list_base[j]/self.scale, torch.cat((LL/self.scale, HL_list[j]/self.scale, LH_list_base[j]/self.scale, self.down[j][1](base_rec)),1))
            sign = torch.sign(LH_list[j])
            tmp = LH_list[j] * sign + 0.5
            tmp[sign == 0] = 0
            LH_list[j] = tmp * sign
            LH_list[j] = (LH_list[j] * self.scale)
            
            # bits = bits + self.coding_HH_list_E[j](HH_list[j]-HH_list_base[j]/self.scale, torch.cat((LL/self.scale, HL_list[j]/self.scale, LH_list[j]/self.scale, HH_list_base[j]/self.scale, self.down[j][2](base_rec)),1))
            sign = torch.sign(HH_list[j])
            tmp = HH_list[j] * sign + 0.5
            tmp[sign == 0] = 0
            HH_list[j] = tmp * sign
            HH_list[j] = (HH_list[j] * self.scale)

            LL = self.wavelet_transform.inverse_trans(LL, HL_list[j], LH_list[j], HH_list[j])
        
        if not train:
            LL = LL[:, :, :height, :width]
        bits = torch.cuda.FloatTensor((1,1,1,1))
        return LL, bits, self.scale*torch.ones(1, device="cuda")
