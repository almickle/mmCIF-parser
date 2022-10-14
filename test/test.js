import { readFileSync } from 'fs'
import 'util'
import { inspect } from 'util'

const PDB = readFileSync('./assets/1axc.cif', 'utf-8').split('\n')
const ATOMS = PDB.filter((line) => line.startsWith('ATOM')).map(array => array.split(' ').filter(element => element !== ''))
const atomLabels = PDB.filter((line) => line.startsWith('_atom_site')).filter((line) => line.startsWith('_atom_sites.fract_transf') === false && line.startsWith('_atom_sites.entry_id') === false)
const customLabels = ['', 'id', 'atom', 'atom_type', '', 'residue', 'chain', 'entity_index', 'residue_index', '', 'x', 'y', 'z', 'occupancy', 'isotropic_temperature_factor', 'formal_charge', 'author_residue_index', 'author_residue', 'author_chain', 'author_atom_type', '' ]
const test = ATOMS[0].map((entry, index) => {
    const label = customLabels[index]
    if(label !== '')
    return (
        {[label]: entry}
    )}).filter((entry) => entry)

    const newAtoms = ATOMS.map((array) => {
        let object = {}
        let author_entries = {}
        array.forEach((entry, index) => {
            const label = customLabels[index]
            if(label !== '') {
                if(label.includes('author')) {
                    if( label.includes('index')) {
                        const int = parseInt(entry)
                        author_entries = ({...author_entries, [label.replace('author_', '')]: int})
                    }
                    else {
                        author_entries = {...author_entries, [label.replace('author_', '')]: entry}
                    }
                    object = {...object, author_entries}
                    return
                }
                if(label === 'id' || label.includes('index') || label === 'formal_charge') {
                    const int = parseInt(entry)
                    object = ({...object, [label]: int})
                    return
                }
                if(label === 'isotropic_temperature_factor' || label === 'x' || label === 'y' || label === 'z' || label === 'occupancy') {
                    const float = parseFloat(entry)
                    object = ({...object, [label]: float})
                    return
                }
                else {
                    object = ({...object, [label]: entry})
                }
            }
        })
        return (
            object
        )
    })

    console.log(newAtoms[0])
    