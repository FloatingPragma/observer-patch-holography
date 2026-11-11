"""Keep explicitly frozen test identities in their authenticated source domain."""


def pytest_configure(config):
    from model.retained_identity_execution.plugin import install
    install(config)
    from model.documentary_release_successor.plugin import install as install_documentary
    install_documentary(config)
