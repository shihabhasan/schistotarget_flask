import sys, os, subprocess, csv

def construct_line( label, line ):
	new_line = []
	if float( label ) == 0.0:
		label = "0"
	new_line.append( label )
	
	for i, item in enumerate( line ):
		if float( item ) == 0.0:
			continue	# sparse!!!
		new_item = "%s:%s" % ( i + 1, item )
		new_line.append( new_item )
	new_line = " ".join( new_line )
	new_line += "\n"
	return new_line

def csv2libsvm(parameter_file):
    input_file=open(parameter_file, "r")
    libsvm_file= open(parameter_file+".libsvm", "w")
    reader = csv.reader(input_file)
    for line in reader:
        new_line = construct_line( 0, line )
        libsvm_file.write( new_line )
    input_file.close()
    libsvm_file.close()

def libsvm_surface(parameter_file):
    csv2libsvm(parameter_file)
    svmpredict="~/schistoprot/tools/libsvm-3.18/svm-predict"
    svmscale="~/schistoprot/tools/libsvm-3.18/svm-scale"
    command1=svmscale+" -l 0 -u 1 "+parameter_file+".libsvm > "+parameter_file+".scale"
    return_code1 = subprocess.call(command1, shell=(sys.platform!="Linux"))
    command2=svmpredict+" "+parameter_file+".scale"+" libsvm_surface.scale.model "+parameter_file+".predict"
    return_code2 = subprocess.call(command2, shell=(sys.platform!="Linux"))

def libsvm_secretory(parameter_file):
    csv2libsvm(parameter_file)
    svmpredict="~/schistoprot/tools/libsvm-3.18/svm-predict"
    svmscale="~/schistoprot/tools/libsvm-3.18/svm-scale"
    command1=svmscale+" -l 0 -u 1 "+parameter_file+".libsvm > "+parameter_file+".scale"
    return_code1 = subprocess.call(command1, shell=(sys.platform!="Linux"))
    command2=svmpredict+" "+parameter_file+".scale"+" libsvm_secretory.scale.model "+parameter_file+".predict"
    return_code2 = subprocess.call(command2, shell=(sys.platform!="Linux"))
