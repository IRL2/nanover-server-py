import sys
import tomllib
import tomli_w

if __name__ == "__main__":
    path = sys.argv[1]

    with open(path, "rb") as infile:
        toml = tomllib.load(infile)

    toml["project"]["requires-python"] = toml["project"]["requires-python"].replace(">=", "==")
    toml["project"]["dependencies"] = [dep.replace(">=", "==") for dep in toml["project"]["dependencies"]]

    with open(path, "wb") as outfile:
        tomli_w.dump(toml, outfile)
