#!/bin/bash
# Evaluate all plates (perform image-based profiling audit)
audit_config="audit_config.yaml"
profile_dir="../0.generate-profiles"

python audit.py \
  --config $audit_config \
  --profile_dir $profile_dir"/profiles" \
  --output_dir "results" \
  --figure_dir "figures"

Rscript --vanilla audit-plate-effects.R \
  --config $audit_config \
  --profile_dir $profile_dir
