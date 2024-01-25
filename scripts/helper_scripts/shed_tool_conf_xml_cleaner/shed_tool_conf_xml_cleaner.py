#!/usr/bin/env python
# Example usage: python shed_tool_conf_xml_cleaner.py
#     -i shed_tool_conf.xml
#     -o shed_tool_conf.xml

import xml.etree.ElementTree as ET
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
from xml.dom import minidom


def merge_xml_sections(input_file, output_file):
    section_counter = 0
    initial_section_number = 0
    tool_counter = 0
    tree = ET.parse(input_file)
    root = tree.getroot()
    children = list(root)
    section_list = []
    for child in children:
        initial_section_number += 1
        section_list.append(tuple(sorted(list(child.items()))))
        # unique sections tuple formatted
        section_list = list(set(section_list))
    toolbox = ET.Element("toolbox")
    toolbox.attrib = root.attrib
    for section in section_list:
        ET.SubElement(toolbox, "section", id=section[0][1],
                      name=section[1][1], version=section[2][1])
    for child in toolbox.iter('section'):
        section_counter += 1
        for originchild in root.iter('section'):
            if child.attrib == originchild.attrib:
                for item in originchild.iter("tool"):
                    tool_counter += 1
                    if 'labels' in item.attrib:
                        tool = ET.SubElement(child, "tool",
                                             file=item.attrib['file'],
                                             guid=item.attrib['guid'],
                                             labels=item.attrib['labels'])
                    else:
                        tool = ET.SubElement(child, "tool",
                                             file=item.attrib['file'],
                                             guid=item.attrib['guid'])
                    for subitem in list(item):
                        buffer = ET.SubElement(tool, subitem.tag)
                        buffer.text = subitem.text
    newtree = ET.ElementTree(toolbox)
    root_newtree = newtree.getroot()
    xmlstr = minidom.parseString(ET.tostring(root_newtree)).toprettyxml(
                                 indent="   ")
    with open(output_file, "w") as f:
        f.write(xmlstr)
    return initial_section_number, section_counter, tool_counter


def xmlbackup(infile):
    with open(infile, 'r') as f:
        with open(infile + '.bak', 'w') as out:
            out.write(f.read())
            print("- Original xml file was saved as .bak")


def _parse_cli_options():
    """
    Parse command line options, returning `parse_args` from `ArgumentParser`.
    """
    parser = ArgumentParser(
        formatter_class=RawTextHelpFormatter,
        description='shed_tool_conf.xml configuration file becomes '
                    'complex over time and tool installations.\n'
                    'Thus, <section> with same attributes id and name \n'
                    '(identity of submenus) are sparse '
                    'in  shed_tool_conf.xml,\nwhich complicate the management '
                    'of submenus in the Galaxy tool panel.\n'
                    'This tool defragment the <section> tags.\n'
                    'After cleaning, each submenu in the Galaxy tool panel '
                    'have a single <section>\n in shed_tool_conf.xml whose '
                    '<tool> children correspond to the tools hosted\n'
                    'in this submenu.\n'
                    'Note that the tools also backup the original '
                    'shed_tool_conf.xml before xml manipulation.',
        usage=" python %(prog)s <options>")
    parser.add_argument("-i", "--input",
                        dest="input_xml",
                        required=True,
                        help="shed_tool_conf.xml file to clean")
    parser.add_argument("-o", "--output_xml",
                        required=True,
                        dest="output_xml",
                        help="Output file for cleaned xml")
    return parser.parse_args()


def __main__():
    args = _parse_cli_options()
    input_xml = args.input_xml
    output_xml = args.output_xml
    xmlbackup(input_xml)
    Init_section, Final_sections, tools = merge_xml_sections(input_xml,
                                                             output_xml)
    print("- %s tools in %s sections were re-packed in %s sections" %
          (tools, Init_section, Final_sections))


if __name__ == "__main__":
    __main__()
