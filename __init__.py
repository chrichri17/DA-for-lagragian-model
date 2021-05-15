# import sys
# from pathlib import Path

# print('start')
# # Avoid import error
# file = Path(__file__).resolve()
# parent, root = file.parent, file.parents[1]
# print(parent, root)
# sys.path.append(str(root))
# # Additionally remove the current file's directory from sys.path
# try:
#     sys.path.remove(str(parent))
# except ValueError: # Already removed
#     pass
