import os as _os
user_path = '/'.join([_os.getenv('HOME'),
                      'aospy_user/aospy_user']).replace('//', '/')
del _os
