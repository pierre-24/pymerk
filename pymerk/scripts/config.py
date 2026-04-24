import argparse
import pathlib
import tomllib

from pymerk.scripts import Config


def main():
    parser = argparse.ArgumentParser(description='pymerk configuration')
    parser.add_argument('-i', '--input', type=pathlib.Path, help='TOML config file')

    args = parser.parse_args()

    # Load Config
    config_data = {}
    if args.input:
        with args.input.open('rb') as f:
            config_data = tomllib.load(f)
    config = Config(**config_data)

    # dump it as TOML
    model = config.model_dump()
    for title, content in model.items():
        print(f'[{title}]')
        for var, val in content.items():
            if val is not None and val != '':
                if isinstance(val, bool):
                    print(f'{var} = {str(val).lower()}')
                else:
                    print(f'{var} = {repr(val)}')


if __name__ == '__main__':
    main()
