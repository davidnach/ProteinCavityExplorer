import { Component } from '@angular/core';
import {Router} from '@angular/router';
import {ServerService} from './server.service';
import {UserInputService} from './userinput.service';
@Component({
  selector: 'app-root',
  templateUrl: './app.component.html',
  styleUrls: ['./app.component.css']
})
export class AppComponent {
  constructor(private ServerService: ServerService,
  private router: Router,private userInputService: UserInputService ){
    this.userInputService.pidEntered.subscribe(
      (status : string) => this.router.navigate(['experiment']));
  }
  title = 'Protein Cavity Explorer';
  experiments = [
    {
    pid: '32',
    name: 'david'}
  ];


  onSubmit(){
    this.ServerService.createExperiment({pid: 'lake'})
    .subscribe((response) => {
      console.log(response);
      console.log(response['pid'])},
    (error) => {console.log(error)});
  }

  submitPid(pid : string){
    console.log(pid);
    this.router.navigate(['experiment']);
  }


  submitExperimentId(experimentId : string){
    console.log(experimentId);

  }

}
