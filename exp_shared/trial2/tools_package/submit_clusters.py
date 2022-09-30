import os
import time

USER="twarowski"

def queueQuery(USER):
  output = os.system(f"qstat -u {USER} > current_jobs.txt")
  print(f"queue query exited with a code {output}")

  jobs_file = open("current_jobs.txt", 'r')
  jobs = jobs_file.readlines()
  jobs_file.close()
  return(len(jobs)-2)


output = os.system(f"ls -l clusterArrayT*.sh  > current_arrays.txt")
print(f"arraysearch exited with a code {output}")
arrs_file = open("current_arrays.txt", 'r')
arrs = arrs_file.readlines()
arrs_file.close()

for script in arrs:
  script = script.split()[-1]
  print(f"current array file to submit: {script}")


  while True:
    job_list_length=queueQuery(USER)
    print(f"the length of a job list for {USER} is {job_list_length}")
    if job_list_length < 1:
      output = os.system(f"qsub {script}")      
      break

    else:
      print(f"query again in 3 min")
      time.sleep(60)
      #print(f"query again in 2 min")
      time.sleep(60)
      #print(f"query again in 1 min")
      time.sleep(60)

        
