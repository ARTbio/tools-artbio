from argparse import ArgumentParser

from bioblend import galaxy


def extract_users(user_list, all_users):
    selected_users = []
    for old_user in user_list:
        for user in all_users:
            if user[u'email'] == old_user:
                selected_users.append(user)
                continue
    return selected_users


def get_user_apikeys(galaxy_instance, users):
    selected_api_keys = []
    for user in users:
        api_key = galaxy_instance.users.get_user_apikey(user[u'id'])
        if api_key == 'Not available.':
            api_key = galaxy_instance.users.create_user_apikey(user[u'id'])
        selected_api_keys.append(api_key)
    return selected_api_keys


def export_user_histories(api_key, GalaxyURL):
    gui = galaxy.GalaxyInstance(url=GalaxyURL, key=api_key)
    histories = gui.histories.get_histories()
    for history in histories:
        print("exporting history %s" % history[u'id'])
        print(gui.histories.export_history(history[u'id'], include_hidden=True))


def _parse_cli_options():
    """
    Parse command line options, returning `parse_args` from `ArgumentParser`.
    """
    parser = ArgumentParser(
        description='bioblend-managed export of histories from a provided \
                     list of galaxy user emails',
        usage=" python %(prog)s <options>")
    parser.add_argument("-g", "--admin-galaxy-url",
                        dest="admin_galaxy_url",
                        required=True,
                        help="Galaxy url of an admin")
    parser.add_argument("-a", "--api-key",
                        required=True,
                        dest="api_key",
                        help="Admin API key")
    parser.add_argument("-e", "--emails",
                        required=True,
                        dest="emails",
                        nargs='+',
                        help="A list of user emails")
    return parser.parse_args()


def __main__():
    args = _parse_cli_options()
    ADMIN_KEY = args.api_key
    GalaxyURL = args.admin_galaxy_url
    gi = galaxy.GalaxyInstance(url=GalaxyURL, key=ADMIN_KEY)
    all_users = gi.users.get_users()
    selected_users = extract_users(args.emails, all_users)
    user_api_keys = get_user_apikeys(gi, selected_users)
    for api_key in user_api_keys:
        export_user_histories(api_key, GalaxyURL)


if __name__ == "__main__":
    __main__()

