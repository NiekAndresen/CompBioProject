with open('/home/niek/HSA_data/____data_experiment_1_2', 'r+') as f:
    head = f.readline()
    first = True
    line = head.split(',')
    newline = line.copy()
    for i,word in zip(range(len(line)), line):
        if word == 'OpenModWindow':
            print("found OpenWindow")
            if first:
                newline[i] = 'OpenModWindow1'
                first = False
            else:
                newline[i] = 'OpenModWindow2'
    f.seek(0)
    newhead = ''.join([w + ',' for w in newline])[:-1]
    print(newhead)
    f.write(newhead)
