from click import command, pass_obj


@command("test")
# Pass context object
@pass_obj
def cli():
    return run()


def run():
    pass
