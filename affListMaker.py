import re

file_read = open("./affParamJS.txt", "r")
lines = file_read.readlines()
file_read.close()

file_write = open("./affParam.txt", "w")

for line in lines:
    new_line = re.sub("element", "", line)
    # new_line = new_line.replace(";", ",")
    new_line = new_line.replace("=", "")
    new_line = new_line.replace("'", "")
    new_line = re.sub("([;]).*?([\[])", "\g<1>\g<2>", new_line)
    new_line = re.sub("[\[].*?[\]]", "", new_line)
    new_line = new_line.replace(";", "\t")
    new_line = new_line[:-2]
    file_write.write(new_line + "\n")
    print(new_line)