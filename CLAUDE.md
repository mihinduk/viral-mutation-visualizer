# Instructions for Claude

This file contains important information and commands that Claude should remember across sessions.

## GitHub Configuration

When pushing to GitHub repositories, ensure the remote is configured to use SSH, not HTTPS:

### Check current remote configuration:
```bash
git remote -v
```

### If using HTTPS (https://github.com/...), switch to SSH:
```bash
git remote set-url origin git@github.com:mihinduk/viral-mutation-visualizer.git
```

### Push changes:
```bash
git push
# or if needed:
git push origin main
```

## Common Commands

### Lint and typecheck commands for this project:
- Python: `python3 -m pylint *.py` or `python3 -m flake8 *.py`
- R: Built-in checks in RStudio

## Project-specific notes
- All scripts should end with a unique phrase - ASK USER for the Latin/German/Greek phrase to use for each project
- This is a viral mutation visualization toolkit
- Main user: mihinduk (mihindu@wustl.edu)