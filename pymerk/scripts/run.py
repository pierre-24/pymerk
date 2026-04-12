import argparse
import pathlib
import tomllib

from pymerk.scripts import Config


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=pathlib.Path, help='Input file')

    args = parser.parse_args()

    config = Config()

    if args.input:
        with args.input.open('rb') as f:
            config = Config(**tomllib.load(f))

    print(config.model_dump())


if __name__ == '__main__':
    main()
