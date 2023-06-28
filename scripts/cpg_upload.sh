TOP_LEVEL_FOLDER=s3://imaging-platform/projects/2018_05_30_ResistanceMechanisms_Kapoor
PROJECT_DIRECTORY=cpg0028-kelley-resistance
LIST_OF_BATCHES=batches.txt
PROJECT_NESTING=broad

mkdir -p log # to log the output

# sync the images

parallel \
  -a ${LIST_OF_BATCHES} \
  --max-procs 4 \
  --eta \
  --joblog log/upload.log \
  --results log/upload \
  --files \
  --keep-order \
  aws s3 sync \
  --acl bucket-owner-full-control \
  --metadata-directive REPLACE \
  "${TOP_LEVEL_FOLDER}"/{1}/images/ \
  s3://cellpainting-gallery/${PROJECT_DIRECTORY}/${PROJECT_NESTING}/images/${BATCH}/images/{1}
  
