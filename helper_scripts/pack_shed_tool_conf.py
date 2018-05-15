#!/usr/bin/env python
# Example usage: python pack_shed_tool_conf.py
#     -i shed_tool_conf.xml
#     -o compact_shed_tool_conf.xml

from xml.dom import minidom
import xml.etree.ElementTree as ET
from argparse import ArgumentParser


def merge_xml_sections(input_file, output_file):
    section_counter = 0
    initial_section_number = 0
    tree = ET.parse(input_file)
    root = tree.getroot()
    children = root.getchildren()
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
                    if 'labels' in item.attrib:
                        tool = ET.SubElement(child, "tool",
                                             file=item.attrib['file'],
                                             guid=item.attrib['guid'],
                                             labels=item.attrib['labels'])
                    else:
                        tool = ET.SubElement(child, "tool",
                                             file=item.attrib['file'],
                                             guid=item.attrib['guid'])
                    tool = ET.SubElement(child, "tool",
                                         file=item.attrib['file'],
                                         guid=item.attrib['guid'],
                                         labels=labels)
                    for subitem in item.getchildren():
                        buffer = ET.SubElement(tool, subitem.tag)
                        buffer.text = subitem.text
    newtree = ET.ElementTree(toolbox)
    root_newtree = newtree.getroot()
    xmlstr = minidom.parseString(ET.tostring(root_newtree)).toprettyxml(
                                 indent="   ")
    with open(output_file, "w") as f:
        f.write(xmlstr)
    return initial_section_number, section_counter


def _parse_cli_options():
    """
    Parse command line options, returning `parse_args` from `ArgumentParser`.
    """
    parser = ArgumentParser(usage=" python %(prog)s <options>")
    parser.add_argument("-i", "--input",
                        dest="input_xml",
                        required=True,
                        help="shed_tool_conf.xml or migrated_tool_conf.xml \
                        that needs to be compacted .")
    parser.add_argument("-o", "--output_xml",
                        required=True,
                        dest="output_xml",
                        help="Output file for compacted xml")
    return parser.parse_args()


def __main__():
    args = _parse_cli_options()
    input_xml = args.input_xml
    output_xml = args.output_xml
    Init_section, Final_sections = merge_xml_sections(input_xml, output_xml)
    print("%s section tool_shed_conf.xml was compacted in %s sections" %
          (Init_section, Final_sections))


if __name__ == "__main__":
    __main__()
