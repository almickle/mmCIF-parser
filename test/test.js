import { readFileSync } from 'fs'
import { inspect } from 'util'
import 'util'


const PDB = readFileSync('./assets/1ul1.cif', 'utf-8').split('\n')
const ATOMS = PDB.filter((line) => line.startsWith('ATOM')).map(array => array.split(' ').filter(element => element !== ''))
const atomLabels = PDB.filter((line) => line.startsWith('_atom_site')).filter((line) => line.startsWith('_atom_sites.fract_transf') === false && line.startsWith('_atom_sites.entry_id') === false)
const customLabels = ['', 'id', 'atom', 'atom_type', '', 'residue', 'chain', 'entity_index', 'residue_index', '', 'x', 'y', 'z', 'occupancy', 'isotropic_temperature_factor', 'formal_charge', 'author_residue_index', 'author_residue', 'author_chain', 'author_atom_type', '' ]

let object = {}
    const newAtoms = ATOMS.map((array) => {
        let atom_info = {}
        let author_entries = {}
        array.forEach((entry, index) => {
            const label = customLabels[index]
            if(label !== '') {
                if(label.includes('author')) {
                    if( label.includes('index')) {
                        const int = parseInt(entry)
                        author_entries = {...author_entries, [label.replace('author_', '')]: int}
                    }
                    else {
                        author_entries = {...author_entries, [label.replace('author_', '')]: entry}
                    }
                    atom_info = {...atom_info, author_entries}
                    return
                }
                if(label === 'id' || label.includes('index') || label === 'formal_charge') {
                    const int = parseInt(entry)
                    atom_info = ({...atom_info, [label]: int})
                    return
                }
                if(label === 'isotropic_temperature_factor' || label === 'x' || label === 'y' || label === 'z' || label === 'occupancy') {
                    const float = parseFloat(entry)
                    atom_info = ({...atom_info, [label]: float})
                    return
                }
                else {
                    atom_info = ({...atom_info, [label]: entry})
                }
            }
        })
        return (
            atom_info
        )
    })
    const entityStart = PDB.indexOf('_entity.details ')
    const chainData = []

    for(let i=entityStart; i <= 1000; i++) {
        if(PDB[i] === '# ') {
            break
        } else {
            chainData.push(PDB[i])
        }
    }

    chainData.shift()

    const chainInfo = chainData.map((line) => replaceSpaces(line, 0)).map((line) => replaceQuotes(line, 0)).map((line) => line.split(' ')).map((array) => array.filter((element) => element !== ''))

    function replaceSpaces(line, count) {
        if(count === line.length-1) {
            return (
                line
            )
        }
        const newLine = line.replace(/(?<='\w+)\s/, '_')
        count++
        return replaceSpaces(newLine, count)
    }

    function replaceQuotes(line, count) {
        if(count >= 2) {
            return (
                line
            )
        }
        const newLine = line.replace(/['"`]/, '')
        count++
        return replaceQuotes(newLine, count)
    }


    const chainArray = chainInfo.map((entry) => {
        let chainObj = {}
        entry.forEach((entry, index) => {
            switch (index) {
                case 1:
                    chainObj = {...chainObj, type: entry }
                    break;
                case 3:
                    chainObj = {...chainObj, name: entry }
                    break;
                case 4:
                    chainObj = {...chainObj, molecular_weight: parseFloat(entry) }
                    break;
                case 5:
                    chainObj = {...chainObj, quantity: parseInt(entry) }
                    break;
            
                default:
                    break;
            }
        })
        return chainObj
    })


    object = {...object, chain_info: chainArray}
    
    let chainCount = 0
    object.chain_info.filter((entry) => entry.type === 'polymer').forEach((entry) => chainCount += entry.quantity)

    const chains = []
    for(let i=0; i < chainCount; i++) {
        chains.push([])
    }

    const chainLabels = newAtoms.map((atom) => atom.chain).filter((chainID, index, array) => chainID !== array[index+1])

    const chainAtoms = chains.map((chain, index) => newAtoms.filter((atom) => atom.chain === chainLabels[index]))


    object = {...object, chains: chainAtoms}
    object = {...object, atoms: newAtoms}


    const backbones = object.chains.map((chain) => chain.filter((residue) => residue.atom_type === 'CA' || residue.atom_type === 'C' || residue.atom_type === 'N' ))
    object = {...object, backbones: backbones }

    console.log(object.backbones[0].length, object.backbones[1].length, object.backbones[2].length, object.backbones[3].length, object.backbones[4].length, object.backbones[5].length)
    console.log(object.chains[0].length, object.chains[1].length, object.chains[2].length, object.chains[3].length, object.chains[4].length, object.chains[5].length)


    