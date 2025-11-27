import re
from src.query.synonyms import CONDITION_SYNONYMS, EXPERIMENT_SYNONYMS, ORGANISM_SYNONYMS

class QueryExpander:
    """
    Expands queries using synonym dictionaries.
    """
    
    def __init__(self):
        self.synonyms = {}
        # Combine all dictionaries for unified search
        self.synonyms.update(CONDITION_SYNONYMS)
        self.synonyms.update(EXPERIMENT_SYNONYMS)
        self.synonyms.update(ORGANISM_SYNONYMS)
        
    def expand_term(self, term):
        """
        Searches for a term in the dictionaries and returns the list of synonyms.
        If no exact match is found, searches for substring.
        """
        term_lower = term.lower()
        
        # 1. Exact match on keys
        if term_lower in self.synonyms:
            return self.synonyms[term_lower]
            
        # 2. Reverse search (if the term is one of the values)
        for key, values in self.synonyms.items():
            if term_lower in [v.lower() for v in values]:
                return values
                
        # 3. Partial search (if the key is contained in the term)
        # Ex: "drought stress" contains "drought"
        for key, values in self.synonyms.items():
            if key in term_lower:
                return values
                
        return [term]

    def expand_query_string(self, query_string):
        """
        Takes a natural query string (e.g., "drought in tomato")
        and returns a list of expanded terms to use in the LLM prompt.
        """
        expanded_concepts = []
        query_lower = query_string.lower()
        
        # Search for known conditions using regex for whole words
        for key, values in self.synonyms.items():
            # Check if the key is present as a whole word
            # Or if any of the values is present
            
            # Build pattern: \b(key|val1|val2)\b
            # Escape values in case they have special characters
            all_terms = [key] + values
            # Filter very short terms (less than 2 chars) to avoid noise, unless they are keys
            terms_to_check = [re.escape(t.lower()) for t in all_terms if len(t) > 1]
            
            if not terms_to_check:
                continue
                
            pattern = r'\b(' + '|'.join(terms_to_check) + r')\b'
            
            if re.search(pattern, query_lower):
                # Found a known concept
                # Format as OR string for the LLM
                or_string = " OR ".join([f'"{v}"' for v in values])
                expanded_concepts.append(f"Concept '{key}': ({or_string})")
                
        return expanded_concepts
