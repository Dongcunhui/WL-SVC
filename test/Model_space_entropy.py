import torch
from torch.nn import functional as F
import Quant            
import PixelCNN
import Wavelet_trans
import numpy as np
import downsample as downsample

class Model_space_entropy(torch.nn.Module):
    def __init__(self, bitplane_number, train_step):
        super(Model_space_entropy, self).__init__()
        self.bitplane_number = bitplane_number

        self.trans_steps = train_step
        self.wavelet_transform = Wavelet_trans.Wavelet()
        self.quant = Quant.Quant()
        self.dequant = Quant.DeQuant()

        # for first bitplane
        self.coding_LL = PixelCNN.PixelCNN()
        self.coding_HL_list = torch.nn.ModuleList([PixelCNN.PixelCNN_Context(1) for _i in range(self.trans_steps)])
        self.coding_LH_list = torch.nn.ModuleList([PixelCNN.PixelCNN_Context(2) for _i in range(self.trans_steps)])
        self.coding_HH_list = torch.nn.ModuleList([PixelCNN.PixelCNN_Context(3) for _i in range(self.trans_steps)])
        # for p in self.parameters():
            # p.requires_grad = False

        # for other bitplane
        self.coding_LL_E_list = torch.nn.ModuleList([])
        self.coding_HL_list_E_list = torch.nn.ModuleList([])
        self.coding_LH_list_E_list = torch.nn.ModuleList([])
        self.coding_HH_list_E_list = torch.nn.ModuleList([])
        self.down_list = torch.nn.ModuleList([])    
        for bitorder in range(self.bitplane_number-1): 
            self.coding_LL_E = PixelCNN.PixelCNN_Context(16 + 1)
            self.coding_LL_E_list.append(self.coding_LL_E)
            self.coding_HL_list_E = torch.nn.ModuleList([PixelCNN.PixelCNN_Context(1 + 16 + 1) for _i in range(self.trans_steps)])
            self.coding_HL_list_E_list.append(self.coding_HL_list_E)
            self.coding_LH_list_E = torch.nn.ModuleList([PixelCNN.PixelCNN_Context(2 + 16 + 1) for _i in range(self.trans_steps)])
            self.coding_LH_list_E_list.append(self.coding_LH_list_E)
            self.coding_HH_list_E = torch.nn.ModuleList([PixelCNN.PixelCNN_Context(3 + 16 + 1) for _i in range(self.trans_steps)])
            self.coding_HH_list_E_list.append(self.coding_HH_list_E)
        
            self.down_1 = torch.nn.ModuleList([downsample.down_sample(1, 16) for _i in range(3)])
            self.down_2 = torch.nn.ModuleList([downsample.down_sample(2, 16) for _i in range(3)])
            self.down_3 = torch.nn.ModuleList([downsample.down_sample(3, 16) for _i in range(3)])
            self.down_4 = torch.nn.ModuleList([downsample.down_sample(4, 16) for _i in range(4)])
            self.down_list.append(torch.nn.ModuleList([self.down_1, self.down_2, self.down_3, self.down_4]))
        

    def forward(self, x, part_bitplane, discard_number, rec_x=0):
        self.part_bitplane = part_bitplane
        # forward transform
        size = x.size()
        width = size[3]
        height = size[2]
        pad_h = int(np.ceil(height / 16)) * 16 - height
        pad_w = int(np.ceil(width / 16)) * 16 - width
        paddings = (0, pad_w, 0, pad_h)
        x = F.pad(x, paddings, 'replicate')

        if part_bitplane == 0:
            self.scale = 2 ** (8 - part_bitplane)
            LL = x
            HL_list = []
            LH_list = []
            HH_list = []
            for i in range(self.trans_steps):
                LL, HL, LH, HH = self.wavelet_transform.forward_trans(LL)
                HL_list.append(self.quant(HL, self.scale))
                LH_list.append(self.quant(LH, self.scale))
                HH_list.append(self.quant(HH, self.scale))

            LL = self.quant(LL, self.scale)

            bits = self.coding_LL(LL)
            sign = torch.sign(LL)
            tmp = LL * sign + 0.5
            tmp[sign == 0] = 0
            LL = tmp * sign
            LL = (LL * self.scale)

            for i in range(self.trans_steps):
                j = self.trans_steps - 1 - i
                
                bits = bits + self.coding_HL_list[j](HL_list[j],LL)
                sign = torch.sign(HL_list[j])
                tmp = HL_list[j] * sign + 0.5
                tmp[sign == 0] = 0
                HL_list[j] = tmp * sign
                HL_list[j] = (HL_list[j] * self.scale)
                
                bits = bits + self.coding_LH_list[j](LH_list[j], torch.cat((LL, HL_list[j]),1))
                sign = torch.sign(LH_list[j])
                tmp = LH_list[j] * sign + 0.5
                tmp[sign == 0] = 0
                LH_list[j] = tmp * sign
                LH_list[j] = (LH_list[j] * self.scale)
                
                bits = bits + self.coding_HH_list[j](HH_list[j], torch.cat((LL, HL_list[j], LH_list[j]),1))
                sign = torch.sign(HH_list[j])
                tmp = HH_list[j] * sign + 0.5
                tmp[sign == 0] = 0
                HH_list[j] = tmp * sign
                HH_list[j] = (HH_list[j] * self.scale)

                LL = self.wavelet_transform.inverse_trans(LL, HL_list[j], LH_list[j], HH_list[j])
        else:    
            # Simulate the previous bitplane
            size = rec_x.size()
            width = size[3]
            height = size[2]
            pad_h = int(np.ceil(height / 16)) * 16 - height
            pad_w = int(np.ceil(width / 16)) * 16 - width
            paddings = (0, pad_w, 0, pad_h)
            rec_x = F.pad(rec_x, paddings, 'replicate')
            self.scale = 2 ** (self.bitplane_number - part_bitplane)
            LL = rec_x
            HL_list_base = []
            LH_list_base = []
            HH_list_base = []
            for i in range(self.trans_steps):
                LL, HL, LH, HH = self.wavelet_transform.forward_trans(LL)
                HL_list_base.append(self.quant(HL, self.scale) * self.scale)
                LH_list_base.append(self.quant(LH, self.scale) * self.scale)
                HH_list_base.append(self.quant(HH, self.scale) * self.scale)
            LL_base = self.quant(LL, self.scale) * self.scale
            LL = LL_base
            for i in range(self.trans_steps):
                j = self.trans_steps - 1 - i
                LL = self.wavelet_transform.inverse_trans(LL, HL_list_base[j], LH_list_base[j], HH_list_base[j])
            base_rec = LL
            
            # Current bitplane
            self.scale = 2 ** (8 - part_bitplane)
            LL = x
            HL_list = []
            LH_list = []
            HH_list = []
            for i in range(self.trans_steps):
                LL, HL, LH, HH = self.wavelet_transform.forward_trans(LL)
                HL_list.append(self.quant(HL, self.scale))
                LH_list.append(self.quant(LH, self.scale))
                HH_list.append(self.quant(HH, self.scale))
            LL = self.quant(LL, self.scale)

            bits = 0
            if discard_number < 13:
                bits = self.coding_LL_E_list[part_bitplane-1](LL - LL_base / self.scale, torch.cat((self.down_list[part_bitplane-1][3][3](base_rec), LL_base/self.scale), 1))
                sign = torch.sign(LL)
                tmp = LL * sign + 0.5
                tmp[sign == 0] = 0
                LL = tmp * sign
                LL = (LL * self.scale)
            else:
                LL = torch.round(LL_base / self.scale)
                sign = torch.sign(LL)
                tmp = LL * sign + self.scale
                tmp[sign == 0] = 0
                LL = tmp * sign

            for i in range(self.trans_steps):
                j = self.trans_steps - 1 - i

                if discard_number < j*3+2+1:
                    bits = bits + self.coding_HL_list_E_list[part_bitplane-1][j](HL_list[j]-HL_list_base[j]/self.scale,torch.cat((LL/self.scale, HL_list_base[j]/self.scale, self.down_list[part_bitplane-1][j][0](base_rec)),1))
                    sign = torch.sign(HL_list[j])
                    tmp = HL_list[j] * sign + 0.5
                    tmp[sign == 0] = 0
                    HL_list[j] = tmp * sign
                    HL_list[j] = (HL_list[j] * self.scale)
                else:
                    HL_list[j] = torch.round(HL_list_base[j])
                    sign = torch.sign(HL_list[j])
                    tmp = HL_list[j] * sign + self.scale
                    tmp[sign == 0] = 0
                    HL_list[j] = tmp * sign



                if discard_number < j * 3 + 1 + 1:
                    bits = bits + self.coding_LH_list_E_list[part_bitplane-1][j](LH_list[j]-LH_list_base[j]/self.scale, torch.cat((LL/self.scale, HL_list[j]/self.scale, LH_list_base[j]/self.scale, self.down_list[part_bitplane-1][j][1](base_rec)),1))
                    sign = torch.sign(LH_list[j])
                    tmp = LH_list[j] * sign + 0.5
                    tmp[sign == 0] = 0
                    LH_list[j] = tmp * sign
                    LH_list[j] = (LH_list[j] * self.scale)
                else:
                    LH_list[j] = torch.round(LH_list_base[j])
                    sign = torch.sign(LH_list[j])
                    tmp = LH_list[j] * sign + self.scale
                    tmp[sign == 0] = 0
                    LH_list[j] = tmp * sign

                if discard_number < j * 3 + 0 + 1:
                    bits = bits + self.coding_HH_list_E_list[part_bitplane-1][j](HH_list[j]-HH_list_base[j]/self.scale, torch.cat((LL/self.scale, HL_list[j]/self.scale, LH_list[j]/self.scale, HH_list_base[j]/self.scale, self.down_list[part_bitplane-1][j][2](base_rec)),1))
                    sign = torch.sign(HH_list[j])
                    tmp = HH_list[j] * sign + 0.5
                    tmp[sign == 0] = 0
                    HH_list[j] = tmp * sign
                    HH_list[j] = (HH_list[j] * self.scale)
                else:
                    HH_list[j] = torch.round(HH_list_base[j])
                    sign = torch.sign(HH_list[j])
                    tmp = HH_list[j] * sign + self.scale
                    tmp[sign == 0] = 0
                    HH_list[j] = tmp * sign

                LL = self.wavelet_transform.inverse_trans(LL, HL_list[j], LH_list[j], HH_list[j])

        LL = LL[:, :, :height, :width]
        
        return LL, bits, self.scale*torch.ones(1, device="cuda")
