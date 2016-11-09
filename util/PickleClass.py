from pickle import dump, load


class PickleClass:

    @staticmethod
    def load(file_path: "relative path to file"):
        with open(file_path, "wb") as f:
            return load(f)

    def save(self, file_path: "relative path to safe file"):
        with open(file_path, "wb") as f:
            dump(self, f)
