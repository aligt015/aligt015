#########################
#Exerice 1 Split-tag

event_times = ctab['TIME']
spec_int_dir = output_dir / 'spec_intervals'
spec_int_dir.mkdir(exist_ok=True)
print("Creating and writing split files...")
split_list = [0, 450, 550, max(event_times)]
splittag.splittag(infiles=transit_exp, outroot=f'./output/spec_intervals/{transit_basename}', time_list=split_list)
