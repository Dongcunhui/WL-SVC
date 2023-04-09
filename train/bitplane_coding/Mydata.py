import numpy as np
import os
from torch.utils.data import Dataset
import torch

class MyData(Dataset):
    def __init__(self, input_file, crop_size, train=0):
        super(MyData, self).__init__()
        
        self.train = train
        if self.train:
            self.input_file_path = input_file[:-13]
            self.name_list = []
            file = open(input_file, "r")
            a = file.readline()
            while(len(a)>1):
                self.name_list.append(a[:-1])
                a = file.readline()
        else:
            self.input_file_path = input_file
            name_list = os.listdir(self.input_file_path)
            self.name_list = [name for name in name_list if "frame_mask" in name]
        self.file_len = len(self.name_list)
        self.crop_size = crop_size
        

    def transform(self, img):
        if self.train:
            height_random = np.random.randint(0, img.shape[1] - (self.crop_size - 1))
            width_random = np.random.randint(0, img.shape[2] - (self.crop_size - 1))
            img = np.array(img, dtype=np.float32)[:, height_random:height_random + self.crop_size, width_random:width_random + self.crop_size]
            img = torch.from_numpy(img)
        else:
            img = np.array(img, dtype=np.float32)
            img = torch.from_numpy(img)
        
        return img

    def __getitem__(self, index):
        self.input_mv = np.load(self.input_file_path+self.name_list[index][:-15]+"mv"+self.name_list[index][-5:])
        self.input_frame_mask = np.load(self.input_file_path+self.name_list[index])
        self.input = np.concatenate((self.input_frame_mask, self.input_mv), 0)
        img = self.transform(self.input)
        return img


    def __len__(self):
        return self.file_len

