import torch


class down_sample_block(torch.nn.Module):
    def __init__(self, internal_channel):
        super(down_sample_block, self).__init__()

        self.internal_channel = internal_channel
        self.conv1 = torch.nn.Conv2d(in_channels=self.internal_channel, out_channels=self.internal_channel, kernel_size=3, stride=2, padding=1)
        self.conv2 = torch.nn.Conv2d(in_channels=self.internal_channel, out_channels=self.internal_channel, kernel_size=3, stride=1, padding=1)
        self.conv3 = torch.nn.Conv2d(in_channels=self.internal_channel, out_channels=self.internal_channel, kernel_size=3, stride=1, padding=1)
        self.relu = torch.nn.ReLU()

    def forward(self, x):
        conv1 = self.conv1(x)
        conv2 = self.conv2(self.relu(conv1))
        conv3 = self.conv3(self.relu(conv2))
        return conv3 + conv1

class down_sample(torch.nn.Module):
    def __init__(self, time, internal_channel):
        super(down_sample, self).__init__()

        self.internal_channel = internal_channel
        self.conv0 = torch.nn.Conv2d(in_channels=1, out_channels=self.internal_channel, kernel_size=3, stride=1, padding=1)

        body = [down_sample_block(self.internal_channel) for _i in range(time)]
        self.body = torch.nn.Sequential(*body)

    def forward(self, x):

        conv0 = self.conv0((x))

        out = self.body(conv0)

        return out
