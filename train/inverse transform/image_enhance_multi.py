import torch

class ResBlock(torch.nn.Module):
    def __init__(self, internal_channel=64):
        super(ResBlock, self).__init__()

        self.internal_channel = internal_channel
        self.padding = torch.nn.ReflectionPad2d(1)
        self.conv1 = torch.nn.Conv2d(in_channels=self.internal_channel, out_channels=self.internal_channel, kernel_size=3, stride=1, padding=0)
        self.relu = torch.nn.ReLU()
        self.conv2 = torch.nn.Conv2d(in_channels=self.internal_channel, out_channels=self.internal_channel, kernel_size=3, stride=1, padding=0)


    def forward(self, x):
        out = self.conv1(self.padding(x))
        out = self.relu(out)
        out = self.conv2(self.padding(out))

        return x + out

class PostProcessing(torch.nn.Module):
    def __init__(self, internal_channel=32):
        super(PostProcessing, self).__init__()

        self.internal_channel = internal_channel
        self.padding = torch.nn.ReflectionPad2d(1)
        
        self.conv1_predict = torch.nn.Conv2d(in_channels=1, out_channels=self.internal_channel, kernel_size=3, stride=1, padding=0)
        self.conv1_resi = torch.nn.Conv2d(in_channels=1, out_channels=self.internal_channel, kernel_size=3, stride=1, padding=0)
        body_predict = [ResBlock(self.internal_channel) for _i in range(2)]
        self.body_predict = torch.nn.Sequential(*body_predict)
        body_resi = [ResBlock(self.internal_channel) for _i in range(2)]
        self.body_resi = torch.nn.Sequential(*body_resi)
        
        self.conv1 = torch.nn.Conv2d(in_channels=2*self.internal_channel+1, out_channels=self.internal_channel, kernel_size=3, stride=1, padding=0)
        body = [ResBlock(self.internal_channel) for _i in range(5)]
        self.body = torch.nn.Sequential(*body)

        self.conv2 = torch.nn.Conv2d(in_channels=self.internal_channel, out_channels=self.internal_channel, kernel_size=3, stride=1, padding=0)

        self.conv3 = torch.nn.Conv2d(in_channels=self.internal_channel, out_channels=1, kernel_size=3, stride=1, padding=0)
        
        self.relu = torch.nn.ReLU()
        
    def forward(self, x, predict, resi):
    
        predict_1 = self.relu(self.conv1_predict(self.padding(predict)))
        predict_2 = self.body_predict(predict_1)
        
        resi_1 = self.relu(self.conv1_resi(self.padding(resi)))
        resi_2 = self.body_resi(resi_1)

        conv1 = self.conv1(self.padding(torch.cat((x, predict_2, resi_2),1)))

        out = self.body(conv1)

        out = self.conv2(self.padding(out))

        out = conv1 + out

        out = self.conv3(self.padding(out))

        return x + out
