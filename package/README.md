# mmCIF-parser

This package contains a function useful for the parsing of PBDx/mmCIF files into JSON. Simply import the function 'parse_mmCIF' and pass a text file as an argument. 


import parse_mmCIF from 'mmcif-parser'

const PDB = parse_mmCIF(text)

console.log(PDB)


The function returns the parsed file in the following format:
```json

{
    atoms: [
        {
            id: number, atom: string, atom_type: string, 
            residue: string, residue_index: number, 
            chain: string, entity_index: number, 
            formal_charge: number, isotropic_temperature_factor: number(float), occupancy: number, 
            x: number(float), y: number(float), z: number(float) 
        }
    ],
    chain_info: [
        {
            type: string, 
            name: string, molecular_weight: number(float), quantity: number
        }
    ],
    chains: [
        { atoms of chain }
    ],
    backbones: [
        { atoms of backbone }
    ]
}
```


