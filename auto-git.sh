#!/bin/bash
git add .
git commit -m "Auto-commit on $(date)" || exit 0
git push origin main