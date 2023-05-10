time_stamp_pattern = r'\[\w{3} \w{3} \d{2} \d{2}:\d{2}:\d{2} \d{4}\]'
rule_name_pattern = r'rule (\w+):$'
jobid_pattern = r'jobid: (\d+)'
running_jobs = {}
finished_jobs = {}
# Open output file
with open("output.log", "w") as outfile:
    # iter on output
    for line in iter(proc.stdout.readline, b''):
        # strips line
        line = line.decode('utf-8').rstrip()
        # writes line to log file
        outfile.write(line + "\n")


        # matches line to time stamp and checks if it matches
        timestamp_match = re.search(time_stamp_pattern, line)
        if timestamp_match:
            # A new timestamp means a new job is starting or has finished
            rule_name = None
            job_number = None
            try:
                next_line = next(iter(proc.stdout.readline, b'')).decode('utf-8').rstrip()
                outfile.write(next_line + "\n")
                # try and matches rule name. if true this means new job init
                rule_name_match = re.search(r'rule (\w+):', next_line)
                if rule_name_match:
                    rule_name = rule_name_match.group(1)
                    # so we will search for the job id 
                    # Read lines until we find the job id or reach the next timestamp
                    while True:
                        try:
                            following_line = next(iter(proc.stdout.readline, b'')).decode('utf-8').rstrip()
                            outfile.write(following_line + "\n")
                            if re.search(time_stamp_pattern, following_line):
                                # We reached the next timestamp, stop reading lines
                                break
                            job_id_match = re.search(r'    jobid: (\d+)', following_line)
                            if job_id_match:
                                job_number = job_id_match.group(1)
                                print(f"Job {job_number} of rule {rule_name} is currently running.")
                                running_jobs[job_number] = rule_name
                        except StopIteration:
                            print("No more lines to read from the output.")       
                else:
                    # Check if the line indicates a job has finished
                    finished_match = re.search(r'Finished job (\d+)\.', next_line)
                    if finished_match:
                        finished_job_number = finished_match.group(1)
                        print(f"Finished job Number: {finished_job_number}")
                        finished_jobs[rule_name] = finished_job_number
                        while True:
                            try:
                                following_following_line = next(iter(proc.stdout.readline, b'')).decode('utf-8').rstrip()
                                outfile.write(following_following_line + "\n")
                                if re.search(time_stamp_pattern, following_following_line):
                                    # We reached the next timestamp, stop reading lines
                                    break
                                progress_match = re.search(r'\d+ of \d+ steps \((\d+)%\) done', following_following_line)
                                if progress_match:
                                    progress = progress_match.group(1)
                                    print(f"Job {finished_job_number} of rule {rule_name} has finished with {progress}% progress.")
                            except StopIteration:
                                print("No more lines to read from the output.")
                                break                                            
            except StopIteration:
                print("No more lines to read from the output.")
                break

                # with open("output.log", "w") as outfile:
                #     for line in iter(proc.stdout.readline, b''):
                #         line = line.decode('utf-8').rstrip()
                #         outfile.write(line + "\n")
                #         timestamp_match = re.search(time_stamp_pattern, line)
                #         if timestamp_match:
                #             # A new timestamp means a new job is starting or has finished
                #             rule_name = None
                #             job_number = None
                #             try:
                #                 next_line = next(iter(proc.stdout.readline, b'')).decode('utf-8').rstrip()
                #                 print(f"{YEL} here is the next line{NC} {next_line}")
                #                 outfile.write(next_line + "\n")
                #                 rule_name_match = re.search(r'rule (\w+):', next_line)
                #                 if rule_name_match:
                #                     rule_name = rule_name_match.group(1)
                #                     # Read the next line and check if it matches the job id pattern
                #                     try:
                #                         next_line = next(iter(proc.stdout.readline, b'')).decode('utf-8').rstrip()
                #                         print(f"{YEL} here is the following line{NC} {next_line}")
                #                         job_id_match = re.search(r'    jobid: (\d+)', next_line)
                #                         if job_id_match:
                #                             job_number = job_id_match.group(1)
                #                             print(f"Job {job_number} of rule {rule_name} is currently running.")

                #                     except StopIteration:
                #                         print("No more lines to read from the output.")
                #             except StopIteration:
                #                 print("No more lines to read from the output.")

                #     # Print final message when Snakemake process is done
                #     custom_output = "Snakemake process finished"
                #     print(custom_output)


                # # Open output file
                # with open("output.log", "w") as outfile:
                #     for line in iter(proc.stdout.readline, b''):
                #         line = line.decode('utf-8').rstrip()
                #         outfile.write(line + "\n")
                #         timestamp_match = re.search(time_stamp_pattern, line)
                #         if timestamp_match:
                #             # A new timestamp means a new job is starting or has finished
                #             rule_name = None
                #             job_number = None
                #             try:
                #                 next_line = next(iter(proc.stdout.readline, b'')).decode('utf-8').rstrip()
                #                 outfile.write(next_line + "\n")
                #                 rule_name_match = re.search(r'rule (\w+):', next_line)
                #                 if rule_name_match:
                #                     rule_name = rule_name_match.group(1)
                #                     # Read lines until we find the job id or reach the next timestamp
                #                     while True:
                #                         try:
                #                             next_line = next(iter(proc.stdout.readline, b'')).decode('utf-8').rstrip()
                #                             outfile.write(next_line + "\n")
                #                             if re.search(time_stamp_pattern, next_line):
                #                                 # We reached the next timestamp, stop reading lines
                #                                 break
                #                             job_id_match = re.search(r'    jobid: (\d+)', next_line)
                #                             if job_id_match:
                #                                 job_number = job_id_match.group(1)
                #                                 print(f"Job {job_number} of rule {rule_name} is currently running.")
                #                                 running_jobs[job_number] = rule_name
                #                         except StopIteration:
                #                             print("No more lines to read from the output.")
                #                             break
                #             except StopIteration:
                #                 print("No more lines to read from the output.")
                            
                #     # Print final message when Snakemake process is done
                #     custom_output = "Snakemake process finished"
                #     print(custom_output)





        
        
