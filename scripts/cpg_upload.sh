TOP_LEVEL_FOLDER=s3://imaging-platform/projects/2018_05_30_ResistanceMechanisms_Kapoor
PROJECT_DIRECTORY=cpg0028-kelley-resistance
LIST_OF_BATCHES=batches.txt
PROJECT_NESTING=broad


# 1. images

folder=images
folder_tag=images
mkdir -p log/upload_${folder_tag} # to log the output

parallel \
  -a ${LIST_OF_BATCHES} \
  --max-procs 4 \
  --eta \
  --joblog log/upload_${folder_tag}.log \
  --results log/upload_${folder_tag} \
  --files \
  --keep-order \
  aws s3 sync \
  --acl bucket-owner-full-control \
  --metadata-directive REPLACE \
  "${TOP_LEVEL_FOLDER}"/{1}/ \
  s3://cellpainting-gallery/${PROJECT_DIRECTORY}/${PROJECT_NESTING}/${folder}/{1}/
  
# 2. workspace/{analysis,backend,load_data_csv}

folder=workspace/analysis
#folder=workspace/backend
#folder=workspace/load_data_csv
folder_tag=analysis
#folder_tag=backend
#folder_tag=load_data_csv
mkdir -p log/upload_${folder_tag} # to log the output

parallel \
  -a ${LIST_OF_BATCHES} \
  --max-procs 4 \
  --eta \
  --joblog log/upload_${folder_tag}.log \
  --results log/upload_${folder_tag} \
  --files \
  --keep-order \
  aws s3 sync \
  --acl bucket-owner-full-control \
  --metadata-directive REPLACE \
  "${TOP_LEVEL_FOLDER}"/{1}/ \
  s3://cellpainting-gallery/${PROJECT_DIRECTORY}/${PROJECT_NESTING}/${folder}/{1}/
  
# 3. workspace/metadata


