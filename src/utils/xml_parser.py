import xml.etree.ElementTree as ET

def parse_sra_xml(exp_xml_string):
    """
    Analiza la cadena ExpXml del resumen SRA y devuelve un diccionario de campos.
    Devuelve None si el análisis falla.
    """
    if not exp_xml_string:
        return None

    wrapped_xml = f"<Root>{exp_xml_string}</Root>"
    
    try:
        root = ET.fromstring(wrapped_xml)
    except ET.ParseError as e:
        print(f"Error analizando XML: {e}")
        return None

    data = {}
    
    # Extraer campos
    data['bioproject'] = root.findtext(".//Bioproject", "")
    data['title'] = root.findtext(".//Title", "")
    
    study_elem = root.find(".//Study")
    data['study_name'] = study_elem.attrib.get("name", "") if study_elem is not None else ""
    
    organism_elem = root.find(".//Organism")
    data['organism'] = organism_elem.attrib.get("ScientificName", "") if organism_elem is not None else ""
    
    data['library_strategy'] = root.findtext(".//Library_descriptor/LIBRARY_STRATEGY", "")
    data['library_source'] = root.findtext(".//Library_descriptor/LIBRARY_SOURCE", "")
    data['library_selection'] = root.findtext(".//Library_descriptor/LIBRARY_SELECTION", "")
    
    data['library_layout'] = (
        "PAIRED"
        if root.find(".//Library_descriptor/LIBRARY_LAYOUT/PAIRED") is not None
        else "SINGLE"
    )
    
    data['tissue'] = root.findtext(".//LIBRARY_NAME", "")
    data['biosample'] = root.findtext(".//Biosample", "")
    
    return data
