# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Modify the ``/ID [<...> <...>]`` line in a PDF file.

This is provided to make generated manuscripts (via ``pdflatex``) be
bitwise identical across runs.
"""

from __future__ import print_function

import argparse


SPLIT_TEXT = b"\n/ID [<"


def do_replace(path, new_id):
    with open(path, "rb") as file_obj:
        contents = file_obj.read()

    # Assert that there is exactly one match.
    pre, post = contents.split(SPLIT_TEXT)
    # ID line expected to be of the form:
    #     /ID [<...> <...>]
    _, post = post.split(b"\n", 1)
    # which would just leave `<...> <...>]` after the split.
    new_id_line = u"<{}> <{}>]\n".format(new_id, new_id)
    new_id_line = new_id_line.encode("ascii")

    new_contents = pre + new_id_line + post
    with open(path, "wb") as file_obj:
        file_obj.write(new_contents)

    print("Updated {}".format(path))


def main():
    description = (
        "Modify the `/ID` property in a PDF generated by `pdflatex`."
    )
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--base", required=True, help="Base path for PDF file to be modified."
    )
    parser.add_argument(
        "--id", dest="id_", required=True, help="The prescribed ID to add."
    )

    args = parser.parse_args()
    filename = "{}.pdf".format(args.base)
    do_replace(filename, args.id_)


if __name__ == "__main__":
    main()
