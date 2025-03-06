old_file = "/home/ktt/CLionProjects/PSU-CIDD-Malaria-Simulation-Dev_GPU/input/rwa/rwa_init_pop.asc"
new_file = "/home/ktt/CLionProjects/PSU-CIDD-Malaria-Simulation-Dev_GPU/input/rwa/rwa_init_pop_resized.asc"

with open(old_file) as old, open(new_file, 'w') as new:
    line_count = 0
    for line in old:
        if line_count > 5:
            line = line.strip()
            # print(line)
            row_data = line.split(" ")
            # print(row_data)
            new_row_data = ''
            if row_data != ['']:
                for loc_data in row_data:
                    if ('-' in loc_data) or (loc_data ==  ''):
                        loc_data2 = loc_data
                    else:
                        loc_data2 = int(int(loc_data)*0.1)
                    new_row_data += str(loc_data2)+' '
                new_row_data += '\n'
            new.write(new_row_data)
        else:
            new.write(line)
        line_count += 1

