from click import command


@command("test")
def cli():
    return run()


def run():
    pass
