import sys
import os
import subprocess
from PIL import Image
import glob

def extractImageData():
    distal = []
    proximal = []
    wholeRegion = []
    images = []
    x1 = 10
    y1 = 66
    x2 = 25
    y2 = 230
    for file in glob.glob('DATA/*.jpg'):
        im = Image.open(file)
        #im = im.crop((x1, x2, y1, y2))
        images.append(im)
    
    testIm = images[150]
    print(len(images))
    images[25].show()
    testIm.show()
    

def runSim(parameterFile):
    os.system("mkdir DATA")
    os.system(("./php_simulator_3s " + str(parameterFile)))
    os.system("./php_mkfig_3s DATA/*.bma")
    return

# def fitCurve(imageData):
#     fit = []


# def fitHistogram():
#     lengths = []

# def saveData(baseFileName):

if __name__ == "__main__":
    #runSim("SimNumber_114087.par")
    extractImageData()
