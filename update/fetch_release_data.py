#!/usr/bin/env python
"""Description:
    Fetch release data from Github and output the result in JSON format
"""
import os
import sys
import json
import argparse
from packaging import version
from github import Github

# define data for externalParentId
CITLINK_1 = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8325244/#:~:text=In%20this%20pipeline%2C%20the%20open%2Dcurated%20generic%20human%20GEM%2C%20Human1%20(11)%2C%20was%20used%20as%20a%20template'
EXT_PARENT_ID = {
        'Mouse-GEM-1.0.0': [
            {
                'id': 'Human-GEM-1.5.0',
                'citLink': CITLINK_1
                }
            ],
        'Rat-GEM-1.0.0': [
            {
                'id': 'Human-GEM-1.5.0',
                'citLink': CITLINK_1
                }
            ],
        'Zebrafish-GEM-1.0.0': [
            {
                'id': 'Human-GEM-1.5.0',
                'citLink': CITLINK_1
                }
            ],
        'Fruitfly-GEM-1.0.0': [
            {
                'id': 'Human-GEM-1.5.0',
                'citLink': CITLINK_1
                }
            ],
        'Worm-GEM-1.0.0': [
            {
                'id': 'Human-GEM-1.5.0',
                'citLink': CITLINK_1
                }
            ],
        }

# define PMID for each specific release (if there is any)
RELEASE_PMID_DICT = {
        }

# define PMID for each GEM
GEM_PMID_DICT = {
        'Yeast-GEM': '31395883',
        'Human-GEM': '32209698',
        'Mouse-GEM': '34282017',
        'Rat-GEM': '34282017',
        'Zebrafish-GEM': '34282017',
        'Fruitfly-GEM': '34282017',
        'Worm-GEM': '34282017',
        }


def writefile(content, outfile, mode='w', is_flush=False):
    """Write text to file"""
    try:
        fpout = open(outfile, mode)
        fpout.write(content)
        if is_flush:
            fpout.flush()
        fpout.close()
    except IOError:
        print(f"Failed to write to {outfile}", file=sys.stderr)


def get_count(repo_name, api_token):
    """Get count for genes, reactions and metabolites"""
    # To be implemented
    return {}


def get_validation(repo_name, api_token):
    """Get quality score for the model"""
    # To be implemented
    return {}


def get_release_data(repo_name, git_api, is_add_version=False):
    """Get release data from Github"""
    repo = git_api.get_repo(repo_name)
    releases = repo.get_releases()
    gem_name = os.path.basename(repo_name)
    gem_data_list = []
    for rel in releases:
        if rel.prerelease:  # ignore pre releases
            continue
        version = rel.tag_name.lstrip('v')
        li = version.split('-')
        if len(li) > 1:  # ignore alpha or beta releases
            continue
        gem_data = {}
        rel_id = f"{gem_name}-{version}"
        gem_data['id'] = rel_id
        if rel_id in EXT_PARENT_ID:
            gem_data['externalParentId'] = EXT_PARENT_ID[rel_id]
        else:
            gem_data['externalParentId'] = []
        gem_data['releaseDate'] = rel.published_at.strftime("%Y-%m-%d")
        if is_add_version:
            gem_data['version'] = version
        gem_data['releaseLink'] = rel.html_url
        if rel_id in RELEASE_PMID_DICT:
            gem_data['PMID'] = RELEASE_PMID_DICT[rel_id]
        elif gem_name in GEM_PMID_DICT:
            gem_data['PMID'] = GEM_PMID_DICT[gem_name]
        else:
            gem_data['PMID'] = ''
        gem_data['count'] = get_count(repo_name, git_api)
        gem_data['validation'] = get_validation(repo_name, git_api)
        gem_data_list.append(gem_data)

    # sort the releases in the ascending order
    gem_data_list = gem_data_list[::-1]

    return gem_data_list


def get_integrated_gems(repo_name, git_api):
    """Get the list of integrated GEMS from Github"""
    repo = git_api.get_repo(repo_name)
    file_content = repo.get_contents("integrated-models/integratedModels.json")
    data = json.loads(file_content.decoded_content.decode())
    gem_set = {}
    for item in data:
        gem_set[item['short_name']] = {}
        gem_set[item['short_name']]['version'] = item['version']
        gem_set[item['short_name']]['date'] = item['date']
    return gem_set


def generate_data(git_api, basedir, is_dryrun):
    """Fetch data from Github and output the result in JSON format"""
    owner = "SysBioChalmers"
    gem_set = get_integrated_gems("MetabolicAtlas/data-files", git_api)
    gem_name_list = gem_set.keys()
    for gem_name in gem_name_list:
        repo_name = f"{owner}/{gem_name}"
        gem_release_data = get_release_data(repo_name, git_api)
        outfile = os.path.join(basedir, 'integrated-models', gem_name,
                               'gemRepository.json')
        if not is_dryrun:
            writefile(json.dumps(gem_release_data, indent=2), outfile)
            msg = f"Result file output to {outfile}"
        else:
            msg = f"DryRun: result file will be output to {outfile}"
        print(msg)


def show_updatable_models(git_api, basedir):
    """Show models that needs to be updated"""
    owner = "SysBioChalmers"
    gem_set = get_integrated_gems("MetabolicAtlas/data-files", git_api)
    for gem_name in gem_set:
        repo_name = f"{owner}/{gem_name}"
        gem_release_data = get_release_data(repo_name, git_api, is_add_version=True)
        latest_gem = gem_release_data[-1]
        v1 = gem_set[gem_name]['version']
        v2 = latest_gem['version']
        if version.parse(v1) < version.parse(v2):
            print(f"{gem_name} can be updated: {v1} => {v2}")


def main():
    """main procedure"""
    parser = argparse.ArgumentParser(
        description='Generate data for GEMs repository',
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument('-d', '--dry', dest='is_dryrun', default=False,
                        help='do not save output to file',
                        action='store_true')
    parser.add_argument('-s', '--show-updatable-model',
                        dest='is_show_updatable_model', default=False,
                        help='show updatable models',
                        action='store_true')

    args = parser.parse_args()
    is_dryrun = args.is_dryrun
    is_show_updatable_model = args.is_show_updatable_model

    try:
        api_token = os.environ['GH_TOKEN']
    except KeyError:
        api_token = None

    git_api = Github(api_token)
    rundir = os.path.dirname(sys.argv[0])
    basedir = os.path.realpath(os.path.join(rundir, os.pardir))

    if is_show_updatable_model:
        show_updatable_models(git_api, basedir)
    else:
        generate_data(git_api, basedir, is_dryrun)


if __name__ == "__main__":
    main()
