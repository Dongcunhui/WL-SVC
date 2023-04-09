import torch
from torch.autograd import Variable
import math
import sys
from torch.nn import functional as F

FIR14 = [ -0.0110, 0.0452, -0.1437, 0.8950, 0.2777, -0.0812, 0.0233, -0.0053]
FIR12 = [-0.0105, 0.0465, -0.1525, 0.6165, 0.6165, -0.1525, 0.0465, -0.0105]
FIR34 = [-0.0053, 0.0233, -0.0812, 0.2777, 0.8950, -0.1437, 0.0452, -0.0110]

class FIR_kernel(torch.nn.Module):
    def __init__(self, trainable_set=False):
        super(FIR_kernel, self).__init__()
        self.FIR1_x = torch.nn.Conv2d(1, 1, (1, 8), padding=0, bias=False)
        self.FIR1_x.weight = torch.nn.Parameter(torch.Tensor([[[[FIR14[0], FIR14[1], FIR14[2], FIR14[3],
                                        FIR14[4], FIR14[5], FIR14[6], FIR14[7]]]]]), requires_grad=trainable_set)

        self.FIR2_x = torch.nn.Conv2d(1, 1, (1, 8), padding=0, bias=False)
        self.FIR2_x.weight = torch.nn.Parameter(torch.Tensor([[[[FIR12[0], FIR12[1], FIR12[2], FIR12[3],
                                                               FIR12[4], FIR12[5], FIR12[6], FIR12[7]]]]]),
                                               requires_grad=trainable_set)

        self.FIR3_x = torch.nn.Conv2d(1, 1, (1, 8), padding=0, bias=False)
        self.FIR3_x.weight = torch.nn.Parameter(torch.Tensor([[[[FIR34[0], FIR34[1], FIR34[2], FIR34[3],
                                                               FIR34[4], FIR34[5], FIR34[6], FIR34[7]]]]]),
                                               requires_grad=trainable_set)

        self.FIR1_y = torch.nn.Conv2d(1, 1, (8, 1), padding=0, bias=False)
        self.FIR1_y.weight = torch.nn.Parameter(torch.Tensor([[[[FIR14[0]], [FIR14[1]], [FIR14[2]], [FIR14[3]],
                                                                [FIR14[4]], [FIR14[5]], [FIR14[6]], [FIR14[7]]]]]),
                                                requires_grad=trainable_set)

        self.FIR2_y = torch.nn.Conv2d(1, 1, (8, 1), padding=0, bias=False)
        self.FIR2_y.weight = torch.nn.Parameter(torch.Tensor([[[[FIR12[0]], [FIR12[1]], [FIR12[2]], [FIR12[3]],
                                                                [FIR12[4]], [FIR12[5]], [FIR12[6]], [FIR12[7]]]]]),
                                                requires_grad=trainable_set)

        self.FIR3_y = torch.nn.Conv2d(1, 1, (8, 1), padding=0, bias=False)
        self.FIR3_y.weight = torch.nn.Parameter(torch.Tensor([[[[FIR34[0]], [FIR34[1]], [FIR34[2]], [FIR34[3]],
                                                                [FIR34[4]], [FIR34[5]], [FIR34[6]], [FIR34[7]]]]]),
                                                requires_grad=trainable_set)

    def forward(self, frame):
        paddings = (3, 4, 0, 0)
        pad_frame = F.pad(frame, paddings, "replicate")

        frame_x0_y0 = frame
        frame_x1_y0 = self.FIR1_x(pad_frame)
        frame_x2_y0 = self.FIR2_x(pad_frame)
        frame_x3_y0 = self.FIR3_x(pad_frame)

        paddings = (0, 0, 3, 4)
        pad_frame_x0_y0 = F.pad(frame_x0_y0, paddings, "replicate")
        pad_frame_x1_y0 = F.pad(frame_x1_y0, paddings, "replicate")
        pad_frame_x2_y0 = F.pad(frame_x2_y0, paddings, "replicate")
        pad_frame_x3_y0 = F.pad(frame_x3_y0, paddings, "replicate")
        frame_x0_y1 = self.FIR1_y(pad_frame_x0_y0)
        frame_x1_y1 = self.FIR1_y(pad_frame_x1_y0)
        frame_x2_y1 = self.FIR1_y(pad_frame_x2_y0)
        frame_x3_y1 = self.FIR1_y(pad_frame_x3_y0)

        frame_x0_y2 = self.FIR2_y(pad_frame_x0_y0)
        frame_x1_y2 = self.FIR2_y(pad_frame_x1_y0)
        frame_x2_y2 = self.FIR2_y(pad_frame_x2_y0)
        frame_x3_y2 = self.FIR2_y(pad_frame_x3_y0)

        frame_x0_y3 = self.FIR3_y(pad_frame_x0_y0)
        frame_x1_y3 = self.FIR3_y(pad_frame_x1_y0)
        frame_x2_y3 = self.FIR3_y(pad_frame_x2_y0)
        frame_x3_y3 = self.FIR3_y(pad_frame_x3_y0)

        output = torch.cat((frame_x0_y0, frame_x1_y0, frame_x2_y0, frame_x3_y0,
                   frame_x0_y1, frame_x1_y1, frame_x2_y1, frame_x3_y1,
                   frame_x0_y2, frame_x1_y2, frame_x2_y2, frame_x3_y2,
                   frame_x0_y3, frame_x1_y3, frame_x2_y3, frame_x3_y3 ), dim=1)

        return output