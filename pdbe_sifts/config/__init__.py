import os
from collections.abc import Mapping
from pathlib import Path

from omegaconf import OmegaConf


def Config(additionalFile=None):
    default_conf_file = Path(__file__).parent.absolute().joinpath("config.yaml")
    default = OmegaConf.load(default_conf_file)

    base = {}
    base_conf_file = os.getenv("RELEASE_CONFIG_DEFAULT")
    if base_conf_file:
        base = OmegaConf.load(base_conf_file)

    extra_conf = {}
    extra_conf_config_file = os.getenv("RELEASE_CONFIG")
    if extra_conf_config_file and base_conf_file:
        extra_conf = OmegaConf.load(extra_conf_config_file)

    add_on = {}
    if additionalFile:
        add_on = OmegaConf.load(additionalFile)

    return OmegaConf.merge(default, base, extra_conf, add_on)


def print_pair(d, start_key=""):
    for k, v in d.items():
        if isinstance(v, Mapping):
            new_start_key = k if not start_key else "_".join([start_key, k])
            print_pair(v, start_key=new_start_key)

        else:
            k = k if not start_key else "_".join([start_key, k])
            print(f"export {k.upper()}={v}")  # noqa: T201

def print_config():
    conf = Config()
    conf2 = OmegaConf.to_container(conf, resolve=True)
    print_pair(conf2)

if __name__ == "__main__":
    print_config()
