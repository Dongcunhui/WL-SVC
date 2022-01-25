import torch
import torch.optim as optim
from torch.autograd import Variable
import math
import sys
from torch.nn import functional as F
import numpy as np

class RoundNoGradient(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x):

        sign = torch.sign(x)
        return torch.floor(x*sign)*sign
    @staticmethod
    def backward(ctx, g):

        return g
# Usage: xq1 = RoundNoGradient.apply(x1)

class Quant(torch.nn.Module):
    def __init__(self):
        super(Quant, self).__init__()
    def forward(self, x, scale):
        return RoundNoGradient.apply(x/scale)

class DeQuant(torch.nn.Module):
    def __init__(self):
        super(DeQuant, self).__init__()
    def forward(self, x, scale):
        return x*scale


# class S_NoGradient(torch.autograd.Function):
    # @staticmethod
    # def forward(ctx, x, alpha, eps=1e-3):
        # # ctx.save_for_backward(x, alpha, eps)
        # alpha_bounded = max(alpha, eps)
        # m = torch.floor(x) + 0.5
        # r = x - m
        # z = torch.tanh(alpha_bounded / 2.) * 2.
        # y = m + torch.tanh(alpha_bounded * r) / z
        # if alpha < eps:
            # return x
        # else:
            # return y

    # @staticmethod
    # def backward(ctx, g):
        
        # return g,None, None
        
class soft_round_Quant(torch.nn.Module):
    def __init__(self):
        super(soft_round_Quant, self).__init__()
    
    # def add_noise(self, x):
        # noise = np.random.uniform(-0.5, 0.5, x.size())
        # noise = torch.Tensor(noise).cuda()
        # return x + noise
    def add_noise(self, x):
        shape=x.size()
        noise = torch.cuda.FloatTensor(shape)
        torch.rand(shape, out=noise)
        return x + noise-0.5
    
    def s_a(self, x, alpha):
        # alpha_bounded = max(alpha, eps)
        if alpha < 1e-3:
            print("alpha is too small!")
            exit()
        alpha_bounded =alpha
        m = torch.floor(x) + 0.5
        r = x - m
        z = torch.tanh(alpha_bounded / 2.) * 2.
        y = m + torch.tanh(alpha_bounded * r) / z
        # if alpha < eps:
            # return x
        # else:
        return y
        
        
    def forward(self, x, alpha):
        return self.s_a(self.add_noise(self.s_a(x, alpha)), alpha)
        # return self.Sa(self.add_noise(self.Sa(x, alpha)), alpha)
