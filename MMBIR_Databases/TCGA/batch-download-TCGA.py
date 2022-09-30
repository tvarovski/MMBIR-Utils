import requests
import json
import sys
import pandas as pd
import re
from math import ceil
import os

#PARAMS
project = "TCGA-LUSC"
strategy = "WXS"
format = "BAM"
slice_size = 200
#END PARAMS


def retreiveProjectData(project, strategy, format, output_name, manifest=False):

  fields = [
      "file_name",
      "cases.case_id",
      "cases.samples.sample_type",
      "cases.disease_type",
      "cases.project.project_id",
      "cases.project.project_name",
      "cases.samples.tumor_descriptor",
      "cases.diagnoses.age_at_diagnosis",
      "cases.demographic.vital_status",
      "cases.diagnoses.days_to_last_followup",
      "cases.demographic.days_to_birth",
      "cases.demographic.days_to_death",
      "cases.diagnoses.tumor_grade",
      "cases.diagnoses.tumor_stage",
      #days to sample collection
      "cases.samples.days_to_collection",
      "cases.samples.days_to_sample_procurment",
      "cases.samples.portions.analytes.aliquots.concentration",
      "cases.samples.portions.analytes.concentration",
      "cases.diagnoses.ajcc_pathologic_stage",
      "cases.diagnoses.ajcc_pathologic_n",
      "cases.diagnoses.ajcc_pathologic_m",
      "cases.diagnoses.ajcc_pathologic_t",
      "cases.diagnoses.figo_stage"
      ]

  fields = ",".join(fields)

  files_endpt = "https://api.gdc.cancer.gov/files" #files

  # This set of filters is nested under an 'and' operator.
  filters = {
      "op": "and",
      "content":[
          {
          "op": "in",
          "content":{
              "field": "cases.project.project_id",
              "value": [project]
              }
          },
          {
          "op": "in",
          "content":{
              "field": "files.experimental_strategy",
              "value": [strategy]
              }
          },
          {
          "op": "in",
          "content":{
              "field": "files.data_format",
              "value": [format]
              }
          }
      ]
  }

  # A POST is used, so the filter parameters can be passed directly as a Dict object.
  params = {
      "filters": filters,
      "fields": fields,
      "format": "TSV",
      "size": "10000"
      }

  if manifest:
      params = {
      "filters": filters,
      "fields": fields,
      "format": "TSV",
      "size": "10000",
      "return_type":"manifest"
      }

  # The parameters are passed to 'json' rather than 'params' in this case
  response = requests.post(files_endpt, headers = {"Content-Type": "application/json"}, json = params)
  original_stdout = sys.stdout

  with open(output_name, 'w') as f:
      sys.stdout = f # Change the standard output to the file we created.
      print(response.content.decode("utf-8"))
      sys.stdout = original_stdout # Reset the standard output to its original value

def createManifestSlices(samples_file_metadata, samples_file_manifest, slice_size, output_name_root):

  #open the file with IDs load to pandas df
  df_metadata = pd.read_csv(samples_file_metadata, sep="\t")
  df_metadata = df_metadata.sort_values(by=["cases.0.case_id"])
  id_list = df_metadata["id"].tolist()

  id_list_length = len(id_list)
  slices = ceil(id_list_length/slice_size)

  print(f"Found {id_list_length} file IDs. Creating {slices} manifest files with {slice_size} IDs each")
  
  df_manifest = pd.read_csv(samples_file_manifest, sep="\t")

  slices_dir="manifest_slices"
  
  if not os.path.exists(slices_dir):
    os.makedirs(slices_dir)

  for i in range(slices):
    #slice the table
    file_uuid_list = id_list[(i*slice_size):((i+1)*slice_size)]

    df_manifest_sliced = df_manifest[df_manifest["id"].isin(file_uuid_list)]

    df_manifest_sliced.to_csv(f"{slices_dir}/{output_name_root}-slice-{i}.txt", index=False, sep="\t")

  print("Finished!")


# retreive project metadata
output_name_meta = f'{project}-{strategy}-{format}-metadata.tsv'
retreiveProjectData(project, strategy, format, output_name_meta, manifest=False)

# retreive project manifest
output_name_manifest = f'{project}-{strategy}-{format}-manifest.tsv'
retreiveProjectData(project, strategy, format, output_name_manifest, manifest=True)

# slice the project manifest by ordered metadata case-ids into chunks of size slice_size
samples_file_manifest = f"{output_name_manifest}"
samples_file_metadata = f"{output_name_meta}"

output_name_root = f"{project}-{strategy}-{format}-manifest"
createManifestSlices(samples_file_metadata, samples_file_manifest, slice_size, output_name_root)
