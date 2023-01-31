import os
import time

USER="twarowski"
max_jobs = 300

def queueQuery(USER):
  output = os.system(f"qstat -u {USER} > current_jobs.txt")
  print(f"queue query exited with a code {output}")

  jobs_file = open("current_jobs.txt", 'r')
  jobs = jobs_file.readlines()
  jobs_file.close()
  curr_jobs=0
  for line in jobs:
    if "cluster" in line:
      curr_jobs+=1
  #return(len(jobs)-2)
  return curr_jobs

if __name__ == "__main__":
  
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
      if job_list_length <= max_jobs:
        output = os.system(f"qsub {script} {script[14:-3]}")      
        break

      else:
        time.sleep(30)
        print(f"query again in 1 min")
        time.sleep(30)
        #print(f"query again in 30s min")
        time.sleep(30)

        
