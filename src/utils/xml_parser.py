import xml.etree.ElementTree as ET

def parse_sra_xml(exp_xml_string, runs_xml_string=None):
    """
    Parses the ExpXml string from SRA summary and returns a dictionary of fields.
    Returns None if parsing fails.
    """
    if not exp_xml_string:
        return None

    wrapped_xml = f"<Root>{exp_xml_string}</Root>"
    
    try:
        root = ET.fromstring(wrapped_xml)
    except ET.ParseError as e:
        print(f"Error parsing XML: {e}")
        return None

    data = {}
    
    # Extract fields
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
    
    data['biosample'] = root.findtext(".//Biosample", "")
    
    # Process Runs if available
    runs = []
    if runs_xml_string:
        try:
            # Runs XML usually comes as <Run .../><Run .../>, so we wrap it
            runs_root = ET.fromstring(f"<Runs>{runs_xml_string}</Runs>")
            for run in runs_root.findall("Run"):
                acc = run.attrib.get("acc")
                if acc:
                    runs.append(acc)
        except ET.ParseError:
            pass # Ignore errors in runs for now
            
    data['run_accessions'] = ",".join(runs)
    
    return data

def parse_sra_full_xml(xml_content):
    """
    Parses the raw XML string from efetch (SRA) using ElementTree.
    Returns a list of dictionaries (one per EXPERIMENT_PACKAGE).
    """
    if not xml_content:
        return []

    try:
        # Wrap in a root element if it's a list of packages not under a single root
        # Usually efetch returns a concatenation of EXPERIMENT_PACKAGE, which is invalid XML without a root.
        # But sometimes it returns an EXPERIMENT_PACKAGE_SET.
        # Let's try to parse directly, if fail, wrap.
        try:
            root = ET.fromstring(xml_content)
        except ET.ParseError:
            # Likely multiple root elements
            root = ET.fromstring(f"<Root>{xml_content}</Root>")
            
    except ET.ParseError as e:
        print(f"Error parsing full SRA XML: {e}")
        return []

    results = []
    
    # Iterate over EXPERIMENT_PACKAGE
    # If wrapped, they are children of Root. If EXPERIMENT_PACKAGE_SET, they are children of that.
    # We can just findall .//EXPERIMENT_PACKAGE
    
    packages = root.findall(".//EXPERIMENT_PACKAGE")
    
    for pkg in packages:
        data = {}
        
        # 1. Study Info
        study = pkg.find(".//STUDY")
        if study is not None:
            descriptor = study.find("DESCRIPTOR")
            if descriptor is not None:
                data['study_title'] = descriptor.findtext("STUDY_TITLE", "")
                data['study_abstract'] = descriptor.findtext("STUDY_ABSTRACT", "")
            
            # External IDs (BioProject)
            data['bioproject'] = ""
            identifiers = study.find("IDENTIFIERS")
            if identifiers is not None:
                for ext_id in identifiers.findall("EXTERNAL_ID"):
                    if ext_id.attrib.get("namespace") == "BioProject":
                        data['bioproject'] = ext_id.text
                        break
        
        # 2. Experiment Info
        experiment = pkg.find(".//EXPERIMENT")
        if experiment is not None:
            design = experiment.find("DESIGN")
            if design is not None:
                data['design_description'] = design.findtext("DESIGN_DESCRIPTION", "")
                
                lib_desc = design.find("LIBRARY_DESCRIPTOR")
                if lib_desc is not None:
                    data['library_strategy'] = lib_desc.findtext("LIBRARY_STRATEGY", "")
                    data['library_source'] = lib_desc.findtext("LIBRARY_SOURCE", "")
                    data['library_selection'] = lib_desc.findtext("LIBRARY_SELECTION", "")
                    
                    layout = lib_desc.find("LIBRARY_LAYOUT")
                    if layout is not None:
                        data['library_layout'] = "PAIRED" if layout.find("PAIRED") is not None else "SINGLE"

            platform = experiment.find("PLATFORM")
            if platform is not None:
                # The first child tag is the platform name (e.g. ILLUMINA)
                if len(platform) > 0:
                    data['platform'] = platform[0].tag
        
        # 3. Organism (from Sample)
        sample = pkg.find(".//SAMPLE")
        if sample is not None:
            sample_name = sample.find("SAMPLE_NAME")
            if sample_name is not None:
                data['organism'] = sample_name.findtext("SCIENTIFIC_NAME", "")
        
        # 4. Runs
        run_set = pkg.find(".//RUN_SET")
        runs = []
        if run_set is not None:
            for run in run_set.findall("RUN"):
                acc = run.attrib.get("accession")
                if acc:
                    runs.append(acc)
        data['run_accessions'] = ",".join(runs)
        
        # Combined text for LLM
        data['full_text'] = f"Title: {data.get('study_title', '')}\nAbstract: {data.get('study_abstract', '')}\nDesign: {data.get('design_description', '')}\nOrganism: {data.get('organism', '')}\nStrategy: {data.get('library_strategy', '')}"
        
        results.append(data)
        
    return results
