import torch
from Model_fir_warp import Model_fir_warp

class ResBlock(torch.nn.Module):
    def __init__(self, internal_channel=64):
        super(ResBlock, self).__init__()

        self.internal_channel = internal_channel
        self.padding = torch.nn.ReflectionPad2d(1)
        self.conv1 = torch.nn.Conv2d(in_channels=self.internal_channel, out_channels=self.internal_channel,
                                     kernel_size=3, stride=1, padding=0)
        self.relu = torch.nn.ReLU(inplace=False)
        self.conv2 = torch.nn.Conv2d(in_channels=self.internal_channel, out_channels=self.internal_channel,
                                     kernel_size=3, stride=1, padding=0)

    def forward(self, x):
        out = self.conv1(self.padding(x))
        out = self.relu(out)
        out = self.conv2(self.padding(out))

        return x + out

class Model_train_predict(torch.nn.Module):
    def __init__(self, ResBlock_num=1):
        super(Model_train_predict, self).__init__()

        self.warp = Model_fir_warp()

        self.padding = torch.nn.ReflectionPad2d(1)
        self.internal_channel = 64
        self.conv1 = torch.nn.Conv2d(in_channels=3, out_channels=self.internal_channel,
                                     kernel_size=3, stride=1, padding=0)

        body = [ResBlock(self.internal_channel) for _i in range(ResBlock_num)]
        self.body = torch.nn.Sequential(*body)

        self.conv2 = torch.nn.Conv2d(in_channels=self.internal_channel, out_channels=self.internal_channel, kernel_size=3,
                                 stride=1, padding=0)

        self.conv3 = torch.nn.Conv2d(in_channels=self.internal_channel, out_channels=1, kernel_size=3,
                                 stride=1, padding=0)
        self.relu = torch.nn.ReLU(inplace=False)
        
    def forward(self, x_l, flo_l, mask_l, x_r, flo_r, mask_r, grid_up, grid_org):

        nearest_L, nearest_R, fir_L, fir_R = self.warp(x_l, flo_l, mask_l, x_r, flo_r, mask_r, grid_up, grid_org)
        
        predict1 = (fir_L + fir_R) / 2.0

        # conv1 = self.relu(self.conv1(self.padding(torch.cat((fir_L, fir_R, predict1), dim=1))))
        
        # res1 = self.body(conv1)

        # conv2 = self.relu(self.conv2(self.padding(res1)))

        # conv3 = self.conv3(self.padding(conv2))
        conv3 = 0

        return predict1, predict1, predict1 + conv3

