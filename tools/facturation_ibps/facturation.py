#!/usr/bin/env python2
# -*- coding: utf-8 -*-


import argparse

import re

import bs4
import html5lib
import PIL
import pandas as pd
import openpyxl


def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument('--input', '-i', action='store', type=str,
                            help="input html code to convert to xlsx")
    the_parser.add_argument('--output', '-o', action='store', type=str,
                            help='xlsx converted file')
    args = the_parser.parse_args()
    return args



def main(input_file, output_file):
    """Script de parsing des fichiers de facturation de l'IBPS"""

    # ouverture fichier input
    with open(input_file, 'r') as file_object:
        facture_html = file_object.read()

    # parsing de la date et de la période de facturation
    date = re.search(ur"Paris le (.*?)</p>", facture_html).group(1)
    periode = re.search(ur"de la prestation (.*?)</p>", facture_html).group(1)

    # parsing des tableaux html avec pandas
    facture_parsed = pd.read_html(
        facture_html,
        thousands='',
        decimal='.',
        flavor='bs4')

    adresse = facture_parsed[0].replace(
        ur'Adresse de l\'appel à facturation : ', ur'', regex=True)

    # supression des symboles € (ça fait planter les calculs dans excel sinon)
    elements = facture_parsed[1].replace(ur'\s*€', ur'', regex=True)

    # conversion des noms de colonnes
    elements_col = elements.iloc[0]
    cout_col = elements_col.str.extract(r'(cout.*)',
                                        expand=False).dropna().iloc[0]
    elements = elements.rename(columns=elements_col).drop(
        elements.index[0])

    misc = facture_parsed[3]

    ref = misc.iloc[:,  # récupération de la référence
                    0].str.extract(r'sur le bon de commande :\s*(.*)$',
                                   expand=False).dropna().iloc[0]

    # ouverture fichier output
    facture_output = openpyxl.load_workbook(
        'template_facture.xlsx', data_only=False, keep_vba=False)
    ws = facture_output.worksheets[0]

    # rajout de l'image de SU qui ne survit pas à la conversion
    img = openpyxl.drawing.image.Image('template_SU.jpg')
    img.anchor = "A1"
    ws.add_image(img)

    # ajout des éléments facturés dans le tableau
    element_row = 23
    for i in range(len(elements)):
        element_row += 1
        ws.cell(row=element_row, column=1, value=elements.iloc[i][u'Objet'])
        ws.cell(
            row=element_row,
            column=2,
            value=elements.iloc[i][u'nombre(s)'])
        ws.cell(
            row=element_row,
            column=4,
            value=elements.iloc[i][cout_col])

    # ajout de l'adresse
    address_row = 7
    for i in range(len(adresse)):
        address_row += 1
        ws.cell(row=address_row, column=3, value=adresse.iloc[i, 0])

    # ajout de la référence/période/date
    ws.cell(row=2, column=3, value=ref)
    ws.cell(row=5, column=5, value=periode)
    ws.cell(row=21, column=5, value=date)

    # export fichier output
    facture_output.save(output_file)
    return


if __name__ == '__main__':
    args = Parser()
    main(args.input, args.output)
