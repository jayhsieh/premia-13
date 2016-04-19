"""Default configuration file
"""

import os

# karrigell_dir is in the namespace where this script is run in k_config.py
root_dir = os.path.join(server_dir,"www")
data_dir = os.path.join(server_dir, "data","www")
cgi_dir = os.path.join(root_dir,"cgi-bin")
cache_dir = os.path.join(data_dir, "cache")

# list of user roles allowed to see directory listings
# if url matches a directory ?
# None means any user. Other values can be 'admin','edit','visit'
allow_directory_listing = [None]

# don't serve files with extension in this list
hide_extensions = [".pdl",".cache"]

# don't serve files with path matching regular expressions in this list
ignore = ["/core.*","/package.*","/conf.*","/data.*"]

# logging file
logging_file = None #os.path.join(karrigell_dir,"logs","access.log")
logging_rotate = "hourly"

# Unicode management
# encode_output : a string = the encoding to use to send data back to the 
# client
output_encoding = None #"utf-8"

# language to translate marked strings to
language = None

# indicates if a complete trace should be printed in case of exception
debug = True

# dictionary of aliases
# if alias["foo"] = "/usr/some_folder" then url /foo/bar.py will
# be resolved to file /usr/some_folder/bar.py
alias = {"dat": r"@PREMIA_DATA_DIR@"}

# use gzip to compress text files ?
gzip = True

# these modules will be available in all scripts namespace
global_modules = [] #[os.path.join(root_dir,"demo","tour","aaa.py")]

# these modules will be executed at startup 
startup_modules = []

# these modules will be executed at shutdown 
shutdown_modules = []

# capture request and response in a database ?
capture = False

# maximum number of sessions
max_sessions = 500
