import os
import json
import matplotlib.pyplot as plt


def get_data_files(data_dir, extension = ".json"):

    res = {}
    for (dirpath, dirnames, filenames) in os.walk(data_dir):
        for filename in filenames:
            if filename.endswith(extension):
                file_path = os.sep.join([dirpath, filename])

                sampleName = filename.split("_")[0]
                planeName = filename.split("_")[1].split(".")[0]
                res[sampleName+"_"+planeName] = file_path
    return res

dataFiles = get_data_files("C:\\Users\\felix\\OneDrive\\Desktop\\2D_output")

allData = {}
for fileName, path in dataFiles.items():
    
    with open(path) as json_file:
        data_dict = json.load(json_file)
        allData[fileName] = data_dict[fileName]
        
#plot 002 fwhm
temp = []
fwhm_002 = []

for fileName, meas in dataFiles.items():
    if fileName.split("_")[1] == "CFS002":
        temp.append(allData[fileName]['temp'])
        fwhm_002.append(allData[fileName]['fwhm'])

temp, fwhm_002 = zip(*sorted(zip(temp, fwhm_002)))

plt.plot(temp, fwhm_002,'o')
plt.xlabel("Temperature (\N{DEGREE SIGN}C)")
plt.ylabel("FWHM ($\AA$)")
plt.title("Co\u2082FeSn (002)")
plt.show()

#plot 004 fwhm
temp = []
fwhm_004 = []

for fileName, meas in dataFiles.items():
    if fileName.split("_")[1] == "CFS004":
        temp.append(allData[fileName]['temp'])
        fwhm_004.append(allData[fileName]['fwhm'])

temp, fwhm_004 = zip(*sorted(zip(temp, fwhm_004)))

plt.plot(temp, fwhm_004,'o')
plt.xlabel("Temperature (\N{DEGREE SIGN}C)")
plt.ylabel("FWHM ($\AA$)")
plt.title("Co\u2082FeSn (004)")
plt.show()

#plot 002 lattice constant
temp = []
a_002 = []

for fileName, meas in dataFiles.items():
    if fileName.split("_")[1] == "CFS002":
        temp.append(allData[fileName]['temp'])
        a_002.append(allData[fileName]['a'])

temp, a_002 = zip(*sorted(zip(temp, a_002)))

plt.plot(temp, a_002,'o')
plt.xlabel("Temperature (\N{DEGREE SIGN}C)")
plt.ylabel("a ($\AA$)")
plt.title("Co\u2082FeSn (002)")
plt.show()

#plot 004 lattice constant
temp = []
a_004 = []

for fileName, meas in dataFiles.items():
    if fileName.split("_")[1] == "CFS004":
        temp.append(allData[fileName]['temp'])
        a_004.append(allData[fileName]['a'])

temp, a_004 = zip(*sorted(zip(temp, a_004)))

plt.plot(temp, a_004,'o')
plt.xlabel("Temperature (\N{DEGREE SIGN}C)")
plt.ylabel("a ($\AA$)")
plt.title("Co\u2082FeSn (004)")
plt.show()

#plot area ratio 004/002
temp = []
ratio = []

for fileName, meas in dataFiles.items():
    if fileName.split("_")[1] == "CFS004":
        sampleName = fileName.split("_")[0]
        temp.append(allData[fileName]['temp'])
        ratio.append(allData[sampleName+"_CFS004"]['area'] 
                     / allData[sampleName+"_CFS002"]['area']) 


temp, ratio = zip(*sorted(zip(temp, ratio)))

plt.plot(temp, ratio,'o')
plt.xlabel("Temperature (\N{DEGREE SIGN}C)")
plt.ylabel("Integrated Area Ratio")
plt.title("Co\u2082FeSn A_(004) / A_(002)")
plt.show()





