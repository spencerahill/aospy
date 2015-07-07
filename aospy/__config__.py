import os as _os

default_colormap = 'RdBu'
user_path = '/'.join([_os.getenv('HOME'),
                      'aospy_user/aospy_user']).replace('//', '/')

del _os
