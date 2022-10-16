# mmCIF-parser

This package contains a function useful for the parsing of PBDx/mmCIF files into JSON. Simply import the function "parse_mmCIF" and pass a text file as an argument. 

```js
import parse_mmCIF from "mmcif-parser"

const PDB = parse_mmCIF(text)

console.log(PDB)
```

The function returns the parsed file in the following format:
```json

{
    "atoms": [
        {
            "id": 1, "atom": "N", "atom_type": "N", 
            "residue": "MET", "residue_index": 1, 
            "chain": "A", "entity_index": 1, 
            "formal_charge": null, "isotropic_temperature_factor": 37.65, "occupancy": 1, 
            "x": 36.886, "y": 53.177, "z": 21.887,
            "author_entries": {
                "residue": "MET", "residue_index": 1, 
                "chain": "A", "atom_type": "N"
            }
        }
    ],
   "chain_info": [
        {
            "type": "polymer", 
            "name": "PCNA", 
            "molecular_weight": 28795.752, 
            "quantity": 3
        }
    ],
    "chains": [
        { "A": ["array of all atoms in chain"] }
    ],
    "backbones": [
        { "A": ["array of N, CA, C atoms comprising the peptide chain backbone"] }
    ],
    "torsion_angles": [
        {
            "A": [ { "residue": "PRO", "phi": -75.456, "psi": -43.580 } ]
        }
    ]
}
```

**The torsion angles are currently inaccurate
