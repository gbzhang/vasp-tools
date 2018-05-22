import os


for i in range(0,6):
	folder = '0' + str(i)
	print(folder)
	os.system('cat ' + folder + '/OUTCAR | grep TOTEN | tail -n 1')
