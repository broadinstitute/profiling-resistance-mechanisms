TOP_LEVEL_FOLDER=s3://imaging-platform/projects/2021_04_26_Production/2021_04_26_Batch1/images
LIST_OF_PLATES_FROM_SAME_BATCH=plate_list.txt # each line in this text file contains a plate name (they must not contain a space; please rename if they do)
PROJECT_DIRECTORY=jump
PROJECT_NESTING=source_4
BATCH=2021_04_26_Batch1

mkdir -p log/${BATCH} # to log the output

parallel \
  -a ${LIST_OF_PLATES_FROM_SAME_BATCH} \
  --max-procs 4 \
  --eta \
  --joblog log/${BATCH}/upload.log \
  --results log/${BATCH}/upload \
  --files \
  --keep-order \
  aws s3 sync \
  --profile jump-cp-role \
  --acl bucket-owner-full-control \
  --metadata-directive REPLACE \
  "${TOP_LEVEL_FOLDER}"/{1} \
  s3://cellpainting-gallery/${PROJECT_DIRECTORY}/${PROJECT_NESTING}/images/${BATCH}/images/{1}
  
