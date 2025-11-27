import xml.etree.ElementTree as ET

def parse_sra_xml(exp_xml_string, runs_xml_string=None):
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
    data['platform'] = root.findtext(".//Platform", "")
    data['total_spots'] = root.find(".//Statistics").attrib.get("total_spots", "") if root.find(".//Statistics") is not None else ""
    data['total_bases'] = root.find(".//Statistics").attrib.get("total_bases", "") if root.find(".//Statistics") is not None else ""
    
    study_elem = root.find(".//Study")
    data['study_name'] = study_elem.attrib.get("name", "") if study_elem is not None else ""
    
    submitter_elem = root.find(".//Submitter")
    data['center_name'] = submitter_elem.attrib.get("center_name", "") if submitter_elem is not None else ""
    
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
    
    # Procesar Runs si estan disponibles
    runs = []
    if runs_xml_string:
        try:
            # Runs XML suele venir como <Run .../><Run .../>, asi que lo envolvemos
            runs_root = ET.fromstring(f"<Runs>{runs_xml_string}</Runs>")
            for run in runs_root.findall("Run"):
                acc = run.attrib.get("acc")
                if acc:
                    runs.append(acc)
        except ET.ParseError:
            pass # Ignoramos errores en runs por ahora
            
    data['run_accessions'] = ",".join(runs)
    
    return data
