#!/usr/bin/env python

## Example usage: python clean_tool_conf.py -i /galaxy/config/shed_tool_conf.xml -o clean_shed_tool_conf.xml

import xml.etree.ElementTree as ET
from os import path
from argparse import ArgumentParser


def check_child(root, children, tooldir, removed_tools = []):
    """
    For each child in children, check if child is tool. If it does not, check
    if child is section. If it is, recurse into section.
    If it has a file attribute, check if the path exists, else remove child from root.
    """
    for child in children:
        if child.tag == "section":
            check_child(root = children, 
                        children = child.getchildren(),
                        tooldir = tooldir, 
                        removed_tools = removed_tools)
        elif child.tag == "tool":
            if path.exists( path.join (tooldir, child.attrib["file"])):
                pass
            else:
                children.remove(child)
                removed_tools.append(child.attrib["file"])
    return removed_tools


def _parse_cli_options():
    """
    Parse command line options, returning `parse_args` from `ArgumentParser`.
    """
    parser = ArgumentParser(usage="usage: python %(prog)s <options>")
    parser.add_argument("-i", "--input",
                        dest="input_xml",
                        required=True,
                        help="shed_tool_conf.xml or migrated_tool_conf.xml \
                        that needs to be cleaned from non-existant entries." )
    parser.add_argument("-o", "--output_xml",
                        required=True,
                        dest="output_xml",
                        help="Output file for cleaned xml")
    return parser.parse_args()

def __main__():
    args = _parse_cli_options()
    input_xml=args.input_xml
    output_xml=args.output_xml
    tree = ET.parse(input_xml)
    root = tree.getroot()
    tooldir = root.attrib["tool_path"]
    children = root.getchildren()
    removed_tools = check_child(root, children, tooldir)
    print "tool xml not found for the follwing tools, removing entries from output xml:"
    for tool in removed_tools:
        print tool
    with open(output_xml, "w") as output:
        output.write(ET.tostring(root))

if __name__ == "__main__":
    __main__()
