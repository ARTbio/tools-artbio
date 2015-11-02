from bioblend.galaxy import GalaxyInstance
from tempfile import NamedTemporaryFile
from bioblend.galaxyclient import ConnectionError
from os.path import basename
import argparse
import datetime
import time
import ftplib
import pysftp
import urlparse
import urllib2
import requests
import json


def parse_args():
    args = argparse.ArgumentParser(description="This script uploads a file to galaxy by ftp or sftp")
    args.add_argument("-fqdn", "--fqdn", required=True, 
                      help="Enter the galaxy server's fully qualified domain name. e.g. usegalaxy.org")
    args.add_argument("-prot", "--protocol", default="http", choices=['http', 'https'],
                     help="Choose Galaxy server's web protocol")
    args.add_argument("-a", "--api_key", required=True, 
                      help="Enter an exisiting users' API key.")
    args.add_argument("-ftp-prot", "--ftp_protocol", default="ftp", choices=["ftp", "sftp"], 
                     help="Upload files by ftp or sftp.")
    args.add_argument("-p", "--port", default=21, type=int, help="Upload Port. Typically 21 for FTP")
    args.add_argument("--http_auth_username", type=str, help="Username for http authentification." )
    args.add_argument("--http_auth_password", type=str, help="Password for http authentification." )
    #return args.parse_args(['-fqdn', 'mississippi.snv.jussieu.fr', '-prot', 'https', '-a', 
    #                        'XXX', '--ftp_protocol', 'ftp', '-p', '21'])
    #return args.parse_args(['-fqdn', 'lbcd41.snv.jussieu.fr', '-prot', 'https', '-a', 
    #                        'XXX', '--ftp_protocol', 'sftp', '-p', '2121',
    #                       '--http_auth_username', "ged", '--http_auth_password', "gcn5pcaf"])
    return args.parse_args()

def inject_auth(f, http_auth_username, http_auth_password):
    def auth(*args, **kwargs):
        kwargs["auth"] = (http_auth_username, http_auth_password)
        return f(*args, **kwargs)
    return auth

def make_post_request(self, url, payload, params=None, files_attached=False, **kwargs):
    """
    Make a POST request using the provided ``url`` and ``payload``.
    The ``payload`` must be a dict that contains the request values.
    The payload dict may contain file handles (in which case the files_attached
    flag must be set to true).
    If the ``params`` are not provided, use ``default_params`` class field.
    If params are provided and the provided dict does not have ``key`` key,
    the default ``self.key`` value will be included in what's passed to
    the server via the request.
    The return value will contain the response body as a JSON object.
    """
    if params is not None and params.get('key', False) is False:
        params['key'] = self.key
    else:
        params = self.default_params

    # Compute data, headers, params arguments for request.post,
    # leveraging the requests-toolbelt library if any files have
    # been attached.
    if files_attached:
        payload.update(params)
        payload = MultipartEncoder(fields=payload)
        headers = self.json_headers.copy()
        headers['Content-Type'] = payload.content_type
        post_params = {}
    else:
        payload = json.dumps(payload)
        headers = self.json_headers
        post_params = params

    r = requests.post(url, data=payload, headers=headers,
                      verify=self.verify, params=post_params, **kwargs)
    if r.status_code == 200:
        return r.json()
    # @see self.body for HTTP response body
    raise ConnectionError("Unexpected response from galaxy: %s" %
                          r.status_code, body=r.text)

def make_delete_request(self, url, payload=None, params=None, **kwargs):
    """
    Make a DELETE request using the provided ``url`` and the optional
    arguments.
    The ``payload`` must be a dict that can be converted into a JSON
    object (via ``json.dumps``)
    If the ``params`` are not provided, use ``default_params`` class field.
    If params are provided and the provided dict does not have ``key`` key,
    the default ``self.key`` value will be included in what's passed to
    the server via the request.
    """
    if params is not None and params.get('key', False) is False:
        params['key'] = self.key
    else:
        params = self.default_params
    r = requests.delete(url, verify=self.verify, data=payload, params=params, **kwargs)
    return r

def make_put_request(self, url, payload=None, params=None, **kwargs):
    """
    Make a PUT request using the provided ``url`` with required payload.
    The ``payload`` must be a dict that can be converted into a JSON
    object (via ``json.dumps``)
    """
    if params is not None and params.get('key', False) is False:
        params['key'] = self.key
    else:
        params = self.default_params
    payload = json.dumps(payload)
    r = requests.put(url, verify=self.verify, data=payload, params=params, **kwargs)
    return r

def create_user(url, api_key):
    """
    Create a new local galaxy user of name test-`date`,
    with username as password.
    """
    gi = GalaxyInstance(url, api_key)
    date=str(datetime.datetime.now().date())
    username = "test-{date}".format(date=date)
    user_email = username + "@test.com"
    password = user_email
    userlist = gi.users.get_users()
    if not user_exists(username, userlist):
        gi.users.create_local_user(username, user_email, password)
        userlist = gi.users.get_users()
    api_key = create_user_api_key(gi, username, userlist)
    return api_key, user_email, password


def user_exists(username, userlist):
    return sum([entry["username"]==username for entry in userlist]) == 1


def create_user_api_key(gi, username, userlist):
    userid = [user['id'] for user in userlist if user['username'] == username][0]
    return gi.users.create_user_apikey(userid)
                            
                            
def ftp_upload_file(fqdn, user_email, password, port, tmpfile):
    """
    Use ftplib to upload to galaxy.
    """
    session = ftplib.FTP()
    session.connect(fqdn, port)
    session.login(user_email, password)
    tmpfile_basename = basename(tmpfile)
    with open(tmpfile) as file:
        session.storbinary('STOR {fn}'.format(fn=tmpfile_basename), file)
    session.quit()

    
def sftp_upload_file(fqdn, user_email, password, port, tmpfile):
    """
    Use pysftp to upload files to galaxy.
    """
    with pysftp.Connection(fqdn, username=user_email, password=password, port=port) as sftp:
        sftp.put(tmpfile)
    

def successfull_upload(url, new_api_key, tmpfile):   
    """
    Connect as user, check if tmpfile exists in uploaded ftp files.
    """
    gi = GalaxyInstance(url, new_api_key)
    files = [ True for ftp_file in gi.ftpfiles.get_ftp_files() if str(ftp_file["path"]) == basename(tmpfile)]
    if files:
        return True
    else:
        return False


def main():
    args = parse_args()
    web_protocol, fqdn, ftp_protocol, port, api_key = \
    args.protocol, args.fqdn, args.ftp_protocol, args.port, args.api_key
    if args.http_auth_username:
        GalaxyInstance.make_delete_request = inject_auth(make_delete_request, args.http_auth_username, args.http_auth_password)
        GalaxyInstance.make_post_request = inject_auth(make_post_request, args.http_auth_username, args.http_auth_password)
        GalaxyInstance.make_put_request = inject_auth(make_put_request, args.http_auth_username, args.http_auth_password)
        GalaxyInstance.make_get_request = inject_auth(GalaxyInstance.make_get_request, args.http_auth_username, args.http_auth_password)
    url = web_protocol + "://" + fqdn
    new_api_key, user_email, password = create_user(url, api_key)
    tmpfile = NamedTemporaryFile()
    tmpfile.write("1\n2\n3\n")
    if ftp_protocol == "ftp":
        ftp_upload_file(fqdn, user_email, password, port, tmpfile.name)
    else:
        sftp_upload_file(fqdn, user_email, password, port, tmpfile.name)
    if not successfull_upload(url, new_api_key, tmpfile.name):
        sys.exit("{ftp_protocol} upload to galaxy server {fqdn} failed.".format(ftp_protocol = ftp_protocol, fqdn = fqdn))

if __name__ == "__main__":
    main()
