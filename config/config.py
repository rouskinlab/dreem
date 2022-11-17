import os, yaml

def input_folder():
  return 'input_folder' #TODO

def input_files(extension=['fastq','fq','fasta','fa']):
  return [f for f in os.listdir(input_folder()) if f.split('.')[-1] in extension]

def output_folder():
  return 'output_folder' #TODO

def temp_folder():
  return 'temp_folder' #TODO

def get_instructions():
  return yaml.safe_load(open('config/config.yml'))

def get_config():
  instructions = get_instructions()
  return {
    'input_folder': input_folder(),
    'input_files': input_files(),
    'output_folder': output_folder(),
    'temp_folder': temp_folder(),
    'demultiplexing':{
      'use': instructions['demultiplexing']['use'],
      'input': input_folder(),
      'temp': temp_folder()+'/demultiplexing',
      'output': output_folder()+'/demultiplexing'
    },
    'alignment':{
      'input': output_folder()+'/demultiplexing' if instructions['demultiplexing']['use'] else input_folder(),
      'temp': temp_folder()+'/alignment',
      'output': output_folder()+'/alignment'
    },
    'vectoring':{
      'input': output_folder()+'/alignment',
      'temp': temp_folder()+'/vectoring',
      'output': output_folder()+'/vectoring'
    },
    'clustering':{
      'use': instructions['clustering']['use'],
      'input': output_folder()+'/vectoring',
      'temp': temp_folder()+'/clustering',
      'output': output_folder()+'/clustering'
    },
    'aggregate':{
      'input': output_folder()+'/clustering' if instructions['clustering']['use'] else output_folder()+'/vectoring',
      'temp': temp_folder()+'/aggregate',
      'output': output_folder()+'/aggregate'
    },
    'post_processing':{
      'use': instructions['post_processing']['use'],
      'input': output_folder()+'/aggregate',
      'temp': temp_folder()+'/post_processing',
      'output': output_folder()+'/post_processing'
    }}


  
